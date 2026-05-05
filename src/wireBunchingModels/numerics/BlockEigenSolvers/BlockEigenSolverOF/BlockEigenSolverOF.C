/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BlockEigenSolverOF.H"
#include "denseMatrixHelperFunctions.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "Pstream.H"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSolverOF, 0);
}


// * * * * * * * * * * * * * Local Helper Functions * * * * * * * * * * * * //

namespace
{
    Foam::vector zeroVector()
    {
        return Foam::vector(0, 0, 0);
    }


    void writeVectorToEigen
    (
        Eigen::Matrix<Foam::scalar, Eigen::Dynamic, 1>& vec,
        const Foam::label start,
        const Foam::vector& value
    )
    {
        for (Foam::label i = 0; i < 3; ++i)
        {
            vec(start + i) = value[i];
        }
    }


    void printVector(const Foam::word& name, const Foam::vector& value)
    {
        Foam::Info << name << Foam::nl;
        for (Foam::label i = 0; i < 3; ++i)
        {
            Foam::Info << "    " << value[i] << Foam::nl;
        }
    }
   void printTensor(const Foam::word& name, const Foam::tensor& value)
   {
       Foam::Info << name << Foam::nl
	       << "    " << value.xx() << " " << value.xy() << " " << value.xz() << Foam::nl
	       << "    " << value.yx() << " " << value.yy() << " " << value.yz() << Foam::nl
	       << "    " << value.zx() << " " << value.zy() << " " << value.zz() << Foam::nl;
   }  


    void writeBlockEigenRigidBodyMotion
    (
        const Foam::Time& runTime,
        const Foam::RigidBodySolution& rigidBodySolution,
        const Foam::vector& velocity,
        const Foam::vector& angularMomentum,
        const Foam::vector& acceleration,
        const Foam::vector& torque
    )
    {
        if (!Foam::Pstream::master())
        {
            return;
        }

        static bool headerWritten = false;
        static Foam::label solveIndex = 0;

        const Foam::word startTimeName =
            runTime.timeName(runTime.startTime().value());

        Foam::fileName historyDir;

        if (Foam::Pstream::parRun())
        {
            historyDir =
                runTime.path()
               /".."
               /"postProcessing"
               /"BlockEigenSolve_rigidBodyMotion"
               /startTimeName;
        }
        else
        {
            historyDir =
                runTime.path()
               /"postProcessing"
               /"BlockEigenSolve_rigidBodyMotion"
               /startTimeName;
        }

        Foam::mkDir(historyDir);

        const Foam::fileName historyFile
        (
            historyDir/"BlockEigenSolve_rigidBodyMotion.dat"
        );

        Foam::OFstream os
        (
            historyFile,
            Foam::IOstreamOption(),
            Foam::IOstreamOption::APPEND
        );

        if (!headerWritten)
        {
            os  << "Time" << " "
                << "timeIndex" << " "
                << "solveIndex" << " "
                << "disp_x" << " "
                << "disp_y" << " "
                << "disp_z" << " "
                << "rotCorr_x" << " "
                << "rotCorr_y" << " "
                << "rotCorr_z" << " "
                << "velocity_x" << " "
                << "velocity_y" << " "
                << "velocity_z" << " "
                << "angularMomentum_x" << " "
                << "angularMomentum_y" << " "
                << "angularMomentum_z" << " "
                << "acceleration_x" << " "
                << "acceleration_y" << " "
                << "acceleration_z" << " "
                << "torque_x" << " "
                << "torque_y" << " "
                << "torque_z"
                << Foam::endl;

            headerWritten = true;
        }

        const Foam::vector& displacement = rigidBodySolution.displacement;
        const Foam::vector& rotationCorrection =
            rigidBodySolution.rotationCorrection;

        os  << runTime.value() << " "
            << runTime.timeIndex() << " "
            << solveIndex++ << " "
            << displacement.x() << " "
            << displacement.y() << " "
            << displacement.z() << " "
            << rotationCorrection.x() << " "
            << rotationCorrection.y() << " "
            << rotationCorrection.z() << " "
            << velocity.x() << " "
            << velocity.y() << " "
            << velocity.z() << " "
            << angularMomentum.x() << " "
            << angularMomentum.y() << " "
            << angularMomentum.z() << " "
            << acceleration.x() << " "
            << acceleration.y() << " "
            << acceleration.z() << " "
            << torque.x() << " "
            << torque.y() << " "
            << torque.z()
            << Foam::endl;

        Foam::Info
            << "BlockEigen rigid-body motion written to "
            << historyFile << Foam::endl;
    }
}


// * * * * * * * * * * * * * * Protected Data Functions * * * * * * * * * * //

void Foam::BlockEigenSolverOF::convertFoamMatrixToEigenMatrix
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei,
    Eigen::SparseMatrix<scalar>& A
)
{
    // Colm- making number of rows larger by 6 rows
    const label nRowsOrig = 6*d.size();
    const label nRows = 6*(d.size() + 1);

    Info << "Number of rows in the original matrix = " << nRowsOrig << endl;
    Info << "Number of rows in matrix = " << nRows << endl;

    // Colm- reserving 6*6 = 36 more spaces for coefficients   
    std::vector<Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(36*(d.size() + l.size() + u.size() + 1));

    // -------------------------------------------------------------------------
    // diagonal
    // -------------------------------------------------------------------------

    label globalRowI = 0;

    forAll(d, cellI)
    {
        const scalarSquareMatrix& curD = d[cellI];

        for (label localRowI = 0; localRowI < 6; ++localRowI)
        {
            for (label localColI = 0; localColI < 6; ++localColI)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalRowI + localRowI,
                        globalRowI + localColI,
                        curD(localRowI, localColI)
                    )
                );
            }
        }

        globalRowI += 6;
    }

    // -------------------------------------------------------------------------
    // rigid-body block: 6x6 identity- Colm
    // -------------------------------------------------------------------------

    const label rbRow = 6*d.size();

    for (label i = 0; i < 6; ++i)
    {
        coefficients.push_back
        (
            Eigen::Triplet<scalar>
            (
                rbRow + i,
                rbRow + i,
                1.0
            )
        );
    }

    // -------------------------------------------------------------------------
    // off-diagonal
    // -------------------------------------------------------------------------

    forAll(u, faceI)
    {
        const scalarSquareMatrix& curU = u[faceI];
        const scalarSquareMatrix& curL = l[faceI];

        const label globalRow = 6*own[faceI];
        const label globalCol = 6*nei[faceI];

        // upper
        for (label localRowI = 0; localRowI < 6; ++localRowI)
        {
            for (label localColI = 0; localColI < 6; ++localColI)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalRow + localRowI,
                        globalCol + localColI,
                        curU(localRowI, localColI)
                    )
                );
            }
        }

        // lower
        for (label localRowI = 0; localRowI < 6; ++localRowI)
        {
            for (label localColI = 0; localColI < 6; ++localColI)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalCol + localRowI,
                        globalRow + localColI,
                        curL(localRowI, localColI)
                    )
                );
            }
        }
    }

    // Insert triplets into the matrix
    //label bnRows = rowPointers.size()-1;
    //nRows = blockSize*bnRows;    
    A.resize(nRows, nRows);
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    // Compressing matrix is meant to help performance    
    A.makeCompressed();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BlockEigenSolverOF::BlockEigenSolverOF
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei
)
:
    d_(d),
    l_(l),
    u_(u),
    own_(own),
    nei_(nei)
{}


// ************************************************************************* //

Foam::scalar Foam::BlockEigenSolverOF::solve
(
    Foam::Field<Foam::scalarRectangularMatrix>& foamX,
    const Foam::Field<Foam::scalarRectangularMatrix>& foamB,
    RigidBodySolution& rigidBodySolution,
    const RigidBodyStepData& rigidBodyData,
    const Time& runTime
)
{
    Eigen::SparseMatrix<scalar> A;
    convertFoamMatrixToEigenMatrix(d_, l_, u_, own_, nei_, A);

    const label nRows = A.rows();
    const label rigidStart = nRows - 6; //Colm- counter for rigidBody

    Info << "Number of rows in matrix called by solver = " << nRows << endl;

    // -------------------------------------------------------------------------
    // Build RHS vector
    // -------------------------------------------------------------------------

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(nRows);
    label index = 0;

    forAll(foamB, cellI)
    {
        b(index++) = foamB[cellI](0,0);
        b(index++) = foamB[cellI](1,0);
        b(index++) = foamB[cellI](2,0);
        b(index++) = foamB[cellI](3,0);
        b(index++) = foamB[cellI](4,0);
        b(index++) = foamB[cellI](5,0);
    }

    // Colm: defining variables for rigid body RHS
    const scalar rbNewmarkGamma = 0.5;
    const scalar rbNewmarkBeta = 0.25;
    const scalar deltaT = runTime.deltaTValue();

    // Colm: Read in variables
    const RigidBodyState& prev = rigidBodyData.previous;
    const RigidBodyState& curr = rigidBodyData.current;

    printVector("Rigid Body Previous Displacement", prev.displacement);
    printTensor("Rigid Body Previous Orientation", prev.orientation);
    printVector("Rigid Body Previous Linear Velocity", prev.velocity);
    printVector("Rigid Body Previous Angular Momentum", prev.angularMomentum);
    printVector("Rigid Body Previous Linear Acceleration", prev.acceleration);
    printVector("Rigid Body Previous Torque", prev.torque);
    printVector("Rigid Body Updated Linear Acceleration", curr.acceleration);
    printVector("Rigid Body Updated Torque", curr.torque);

    Foam::vector rbVelocity = zeroVector();
    Foam::vector rbAngularMomentum = zeroVector();

    // Colm: Newmark-Beta equations
     // Colm: Solve for v and pi here
    // Colm: Not sure what to do with result- following sixDofRigidBodyMotion/FvBeamNewmark structure   
    for (label i = 0; i < 3; ++i)
    {
        rbVelocity[i] =
            prev.velocity[i]
          + deltaT
           *(
                rbNewmarkGamma*curr.acceleration[i]
              + (1.0 - rbNewmarkGamma)*prev.acceleration[i]
            );

        rbAngularMomentum[i] =
            prev.angularMomentum[i]
          + deltaT
           *(
                rbNewmarkGamma*curr.torque[i]
              + (1.0 - rbNewmarkGamma)*prev.torque[i]
            );
    }

    printVector("Rigid Body Updated Velocity", rbVelocity);
    printVector("Rigid Body Updated Angular Momentum", rbAngularMomentum);

    // translational part of rigid-body RHS
    for (label i = 0; i < 3; ++i)
    {
        b(rigidStart + i) =
            prev.displacement[i]
          + deltaT*prev.velocity[i]
          + deltaT*deltaT
           *(
                rbNewmarkBeta*curr.acceleration[i]
              + (0.5 - rbNewmarkBeta)*prev.acceleration[i]
            );
    }

    // rotational correction part of rigid-body RHS
    for (label i = 0; i < 3; ++i)
    {
        b(rigidStart + 3 + i) =
            deltaT*prev.angularMomentum[i]
          + deltaT*deltaT
           *(
                rbNewmarkBeta*curr.torque[i]
              + (0.5 - rbNewmarkBeta)*prev.torque[i]
            );
    }

    // -------------------------------------------------------------------------
    // Build initial guess
    // -------------------------------------------------------------------------

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(nRows);
    index = 0;

    forAll(foamX, cellI)
    {
        x(index++) = foamX[cellI](0,0);
        x(index++) = foamX[cellI](1,0);
        x(index++) = foamX[cellI](2,0);
        x(index++) = foamX[cellI](3,0);
        x(index++) = foamX[cellI](4,0);
        x(index++) = foamX[cellI](5,0);
    }

    // Colm: Rigid body solution vector
    for (label i = 0; i < 6; ++i)
    {
        x(rigidStart + i) = 0.0;
    }

    // -------------------------------------------------------------------------
    // Calculate initial residual
    // -------------------------------------------------------------------------

    // Colm: +1 for extra 6 x 6 block for rigid body
    const label nCells = d_.size() + 1;
    scalar initialResidual = 0.0;

    {
        Eigen::Matrix<scalar, Eigen::Dynamic, 1> Ax(nRows);
        Ax = A*x;

        Field<scalarRectangularMatrix> foamAx
        (
            nCells,
            scalarRectangularMatrix(6, 1, 0.0)
        );

        Field<scalarRectangularMatrix> foamRhs
        (
            nCells,
            scalarRectangularMatrix(6, 1, 0.0)
        );

        label k = 0;
        for (label i = 0; i < nCells; ++i)
        {
            for (label j = 0; j < 6; ++j)
            {
                foamAx[i](j,0) = Ax[k];
                foamRhs[i](j,0) = b[k];
                ++k;
            }
        }

        Field<scalarRectangularMatrix> blockR(foamRhs - foamAx);
        initialResidual = sqrt(sum(magSqr(blockR)));
    }

    // -------------------------------------------------------------------------
    // Solve
    // -------------------------------------------------------------------------

    typedef enum
    {
        SparseLU,
        BiCGSTAB,
        GMRES,
        DGMRES,
        MINRES
    } solvers;

    solvers sol = SparseLU;

    switch (sol)
    {
        case SparseLU:
        {
            Eigen::SparseLU
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::COLAMDOrdering<int>
            > solver(A);

            x = solver.solve(b);
            break;
        }

        case BiCGSTAB:
        {
            Eigen::BiCGSTAB
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::IncompleteLUT<scalar>
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        case GMRES:
        {
            Eigen::GMRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::IncompleteLUT<scalar>
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        case DGMRES:
        {
            Eigen::DGMRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::IncompleteLUT<scalar>
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        case MINRES:
        {
            Eigen::MINRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::Lower|Eigen::Upper,
                Eigen::IdentityPreconditioner
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        default:
            FatalErrorIn
            (
                "Foam::scalar Foam::BlockEigenSolverOF::solve(...)"
            )   << "Undefined Eigen solver."
                << abort(FatalError);
    }

    // -------------------------------------------------------------------------
    // Copy solved beam unknowns back to foamX
    // -------------------------------------------------------------------------

    index = 0;
    forAll(foamX, cellI)
    {
        foamX[cellI](0,0) = x(index++);
        foamX[cellI](1,0) = x(index++);
        foamX[cellI](2,0) = x(index++);
        foamX[cellI](3,0) = x(index++);
        foamX[cellI](4,0) = x(index++);
        foamX[cellI](5,0) = x(index++);
    }

    // -------------------------------------------------------------------------
    // Copy rigid-body result into named output
    // -------------------------------------------------------------------------
    // Colm- Codex idea
    
    rigidBodySolution.displacement = zeroVector();
    rigidBodySolution.rotationCorrection = zeroVector();

    for (label i = 0; i < 3; ++i)
    {
        rigidBodySolution.displacement[i] = x(rigidStart + i);
        rigidBodySolution.rotationCorrection[i] = x(rigidStart + 3 + i);
    }

    Info << "Output rigid-body solution:" << endl;
    Info << "___________________________" << endl;
    printVector("Displacement", rigidBodySolution.displacement);
    printVector("Rotation correction", rigidBodySolution.rotationCorrection);
    printVector("Velocity", rbVelocity);
    printVector("Angular momentum", rbAngularMomentum);
    printVector("Acceleration", curr.acceleration);
    printVector("Torque", curr.torque);
    Info << "___________________________" << endl;

    writeBlockEigenRigidBodyMotion
    (
        runTime,
        rigidBodySolution,
        rbVelocity,
        rbAngularMomentum,
        curr.acceleration,
        curr.torque
    );

    return initialResidual;
}
