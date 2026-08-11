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

    }


    void writeBlockEigenRigidBodyForceCoupling
    (
        const Foam::Time& runTime,
        const Foam::RigidBodyForceCoupling& coupling,
        const Foam::vector& accelerationContribution,
        const Foam::vector& torqueContribution
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
               /"BlockEigenSolve_rigidBodyForceCoupling"
               /startTimeName;
        }
        else
        {
            historyDir =
                runTime.path()
               /"postProcessing"
               /"BlockEigenSolve_rigidBodyForceCoupling"
               /startTimeName;
        }

        Foam::mkDir(historyDir);

        const Foam::fileName historyFile
        (
            historyDir/"BlockEigenSolve_rigidBodyForceCoupling.dat"
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
                << "force_x force_y force_z "
                << "moment_x moment_y moment_z "
                << "position_x position_y position_z "
                << "centreOfRotation_x centreOfRotation_y centreOfRotation_z "
                << "accelerationContribution_x "
                << "accelerationContribution_y "
                << "accelerationContribution_z "
                << "torqueContribution_x "
                << "torqueContribution_y "
                << "torqueContribution_z"
                << Foam::endl;

            headerWritten = true;
        }

        os  << runTime.value() << " "
            << runTime.timeIndex() << " "
            << solveIndex++ << " "
            << coupling.force.x() << " "
            << coupling.force.y() << " "
            << coupling.force.z() << " "
            << coupling.moment.x() << " "
            << coupling.moment.y() << " "
            << coupling.moment.z() << " "
            << coupling.position.x() << " "
            << coupling.position.y() << " "
            << coupling.position.z() << " "
            << coupling.centreOfRotation.x() << " "
            << coupling.centreOfRotation.y() << " "
            << coupling.centreOfRotation.z() << " "
            << accelerationContribution.x() << " "
            << accelerationContribution.y() << " "
            << accelerationContribution.z() << " "
            << torqueContribution.x() << " "
            << torqueContribution.y() << " "
            << torqueContribution.z()
            << Foam::endl;

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
    const scalar rigidBodyNewmarkDisplacementScale,
    Eigen::SparseMatrix<scalar>& A
)
{
    // Colm- making number of rows larger by 6 rows
    const label nRows = 6*(d.size() + 1);

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
    // rigid-body block: displacement/rotation update equations
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

    // Off-diagonal coupling blocks between beam-end and rigid-body DOFs.
    if
    (
        rigidBodyKinematicCoupling_
     && rigidBodyAttachmentCell_ >= 0
     && rigidBodyAttachmentCell_ < d.size()
    )
    {
        const label beamRow = 6*rigidBodyAttachmentCell_;
        for (label rowI = 0; rowI < 3; ++rowI)
        {
            for (label colI = 0; colI < 3; ++colI)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        beamRow + rowI,
                        rbRow + colI,
                        rigidBodyTranslationCoeff_(rowI, colI)
                    )
                );

                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        beamRow + rowI,
                        rbRow + 3 + colI,
                        rigidBodyRotationCoeff_(rowI, colI)
                    )
                );
            }
        }

        const bool addReciprocalRigidBodyForceRows =
            rigidBodyForceCoupling_.active;

        if (addReciprocalRigidBodyForceRows && rigidBodyForceCoupling_.active)
        {
            const tensor rbTranslationForceCoeff =
                rigidBodyNewmarkDisplacementScale
               *rigidBodyBeamForceWCoeff_;

            const tensor rbRotationForceCoeff =
                rigidBodyNewmarkDisplacementScale
               *rigidBodyRotationCoeff_;

            const tensor beamDisplacementForceCoeff =
                rigidBodyNewmarkDisplacementScale
               *rigidBodyBeamForceWCoeff_;

            for (label rowI = 0; rowI < 3; ++rowI)
            {
                for (label colI = 0; colI < 3; ++colI)
                {
                    // Newmark displacement response to incremental beam
                    // attachment stretch:
                    //   dt^2*beta/m * K * (dWrbAttach - dWbeamCell)
                    // where K is rigidBodyBeamForceWCoeff_ ~= CQW/L.
                    coefficients.push_back
                    (
                        Eigen::Triplet<scalar>
                        (
                            rbRow + rowI,
                            rbRow + colI,
                            rbTranslationForceCoeff(rowI, colI)
                        )
                    );

                    coefficients.push_back
                    (
                        Eigen::Triplet<scalar>
                        (
                            rbRow + rowI,
                            rbRow + 3 + colI,
                            rbRotationForceCoeff(rowI, colI)
                        )
                    );

                    coefficients.push_back
                    (
                        Eigen::Triplet<scalar>
                        (
                            rbRow + rowI,
                            beamRow + colI,
                           -beamDisplacementForceCoeff(rowI, colI)
                        )
                    );
                }
            }
        }

        Info<< "BlockEigen kinematic coupling terms added: beamRows="
            << beamRow << ".." << beamRow + 2
            << ", rigid translation rows/cols=" << rbRow << ".." << rbRow + 2
            << ", rigid rotation rows/cols=" << rbRow + 3
            << ".." << rbRow + 5
            << ", reciprocal rigid-body rows "
            << (addReciprocalRigidBodyForceRows ? "added" : "inactive")
            << endl;
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
    nei_(nei),
    rigidBodyKinematicCoupling_(false),
    rigidBodyAttachmentCell_(-1),
    rigidBodyTranslationCoeff_(tensor::zero),
    rigidBodyRotationCoeff_(tensor::zero),
    rigidBodyBeamForceWCoeff_(tensor::zero),
    rigidBodyBeamForceThetaCoeff_(tensor::zero),
    rigidBodyMomentArm_(vector::zero),
    rigidBodyAttachmentDisplacementPrevious_(vector::zero),
    rigidBodyForceCoupling_()
{}


Foam::BlockEigenSolverOF::BlockEigenSolverOF
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei,
    const bool rigidBodyKinematicCoupling,
    const label rigidBodyAttachmentCell,
    const tensor& rigidBodyTranslationCoeff,
    const tensor& rigidBodyRotationCoeff,
    const tensor& rigidBodyBeamForceWCoeff,
    const tensor& rigidBodyBeamForceThetaCoeff,
    const vector& rigidBodyMomentArm,
    const vector& rigidBodyAttachmentDisplacementPrevious,
    const RigidBodyForceCoupling& rigidBodyForceCoupling
)
:
    d_(d),
    l_(l),
    u_(u),
    own_(own),
    nei_(nei),
    rigidBodyKinematicCoupling_(rigidBodyKinematicCoupling),
    rigidBodyAttachmentCell_(rigidBodyAttachmentCell),
    rigidBodyTranslationCoeff_(rigidBodyTranslationCoeff),
    rigidBodyRotationCoeff_(rigidBodyRotationCoeff),
    rigidBodyBeamForceWCoeff_(rigidBodyBeamForceWCoeff),
    rigidBodyBeamForceThetaCoeff_(rigidBodyBeamForceThetaCoeff),
    rigidBodyMomentArm_(rigidBodyMomentArm),
    rigidBodyAttachmentDisplacementPrevious_
    (
        rigidBodyAttachmentDisplacementPrevious
    ),
    rigidBodyForceCoupling_(rigidBodyForceCoupling)
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
    // Colm: defining variables for rigid body RHS
    const scalar rbNewmarkGamma = 0.5;
    const scalar rbNewmarkBeta = 0.25;
    const scalar deltaT = runTime.deltaTValue();

    scalar rigidBodyNewmarkDisplacementScale = 0.0;

    if (rigidBodyForceCoupling_.active)
    {
        if (mag(rigidBodyForceCoupling_.mass) < SMALL)
        {
            FatalErrorInFunction
                << "blockEigenForceCoupling requested with zero rigid-body mass"
                << abort(FatalError);
        }

        rigidBodyNewmarkDisplacementScale =
            sqr(deltaT)*rbNewmarkBeta/rigidBodyForceCoupling_.mass;
    }

    Eigen::SparseMatrix<scalar> A;
    convertFoamMatrixToEigenMatrix
    (
        d_,
        l_,
        u_,
        own_,
        nei_,
        rigidBodyNewmarkDisplacementScale,
        A
    );

    const label nRows = A.rows();
    const label rigidStart = nRows - 6; //Colm- counter for rigidBody

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

    // Colm: Read in variables
    const RigidBodyState& prev = rigidBodyData.previous;
    const RigidBodyState& curr = rigidBodyData.current;

    vector currAcceleration = curr.acceleration;
    vector currTorque = curr.torque;

    if (rigidBodyForceCoupling_.active)
    {
        if (mag(rigidBodyForceCoupling_.mass) < SMALL)
        {
            FatalErrorInFunction
                << "blockEigenForceCoupling requested with zero rigid-body mass"
                << abort(FatalError);
        }

        const vector accelerationContribution =
            rigidBodyForceCoupling_.force/rigidBodyForceCoupling_.mass;

        const vector momentArm =
            rigidBodyForceCoupling_.position
          - rigidBodyForceCoupling_.centreOfRotation;

        const vector globalMomentContribution =
            rigidBodyForceCoupling_.moment
          + (momentArm ^ rigidBodyForceCoupling_.force);

        const vector torqueContribution =
            rigidBodyForceCoupling_.orientation.T()
          & globalMomentContribution;

        currAcceleration += accelerationContribution;
        currTorque += torqueContribution;

        writeBlockEigenRigidBodyForceCoupling
        (
            runTime,
            rigidBodyForceCoupling_,
            accelerationContribution,
            torqueContribution
        );
    }

    Foam::vector rbVelocity = Foam::vector::zero;
    Foam::vector rbAngularMomentum = Foam::vector::zero;

    // Colm: Newmark-Beta equations
     // Colm: Solve for v and pi here
    // Colm: Not sure what to do with result- following sixDofRigidBodyMotion/FvBeamNewmark structure   
    for (label i = 0; i < 3; ++i)
    {
        rbVelocity[i] =
            prev.velocity[i]
          + deltaT
           *(
                rbNewmarkGamma*currAcceleration[i]
              + (1.0 - rbNewmarkGamma)*prev.acceleration[i]
            );

        rbAngularMomentum[i] =
            prev.angularMomentum[i]
          + deltaT
           *(
                rbNewmarkGamma*currTorque[i]
              + (1.0 - rbNewmarkGamma)*prev.torque[i]
            );
    }

    // translational part of rigid-body RHS
    for (label i = 0; i < 3; ++i)
    {
        b(rigidStart + i) =
            prev.displacement[i]
          + deltaT*prev.velocity[i]
          + deltaT*deltaT
           *(
                rbNewmarkBeta*currAcceleration[i]
              + (0.5 - rbNewmarkBeta)*prev.acceleration[i]
            )
          - rigidBodyAttachmentDisplacementPrevious_[i];
    }

    // rotational correction part of rigid-body RHS
    for (label i = 0; i < 3; ++i)
    {
        b(rigidStart + 3 + i) =
            deltaT*prev.angularMomentum[i]
          + deltaT*deltaT
           *(
                rbNewmarkBeta*currTorque[i]
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
    
    rigidBodySolution.displacement = Foam::vector::zero;
    rigidBodySolution.rotationCorrection = Foam::vector::zero;
    rigidBodySolution.velocity = rbVelocity;
    rigidBodySolution.angularMomentum = rbAngularMomentum;
    rigidBodySolution.acceleration = currAcceleration;
    rigidBodySolution.torque = currTorque;

    for (label i = 0; i < 3; ++i)
    {
        rigidBodySolution.displacement[i] =
            rigidBodyAttachmentDisplacementPrevious_[i]
          + x(rigidStart + i);
        rigidBodySolution.rotationCorrection[i] = x(rigidStart + 3 + i);
    }

    writeBlockEigenRigidBodyMotion
    (
        runTime,
        rigidBodySolution,
        rbVelocity,
        rbAngularMomentum,
        currAcceleration,
        currTorque
    );

    return initialResidual;
}
