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
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BlockEigenSolverOF.H"

// #include "OFstream.H"
// #include "polyMesh.H"
// #include "addToRunTimeSelectionTable.H"
// #include "fvMesh.H"
// #include "beamModel.H"
// #include "multibeamFvBlockMatrix.H"
#include "denseMatrixHelperFunctions.H"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSolverOF, 0);

    // addToRunTimeSelectionTable
    // (
    //     blockVector6Solver, BlockEigenSolverOF, symMatrix
    // );

    // addToRunTimeSelectionTable
    // (
    //     blockVector6Solver, BlockEigenSolverOF, asymMatrix
    // );
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

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
    // if (BlockLduSolver::debug)
    // {
    //     Info<< this->typeName
    //         << ": copying matrix coefficients into Eigen format"
    //         << endl;
    // }

    // Colm- making number of rows larger by 6 rows
    const label nRows_orig = 6*d.size();
    const label nRows = 6*(d.size() + 1);
    Info << "Number of rows in the original matrix = " << nRows_orig <<endl;
    Info << "Number of rows in matrix = " << nRows << endl;
      
    // Block CSR matrix storage

    //const labelList& rowPointers = mbMatrix.blockRowPointers();

    //const labelList& columnIndices = mbMatrix.blockColumnIndices();
    //const scalarField& coeffs = mbMatrix.blockCoeffs();

    //label blockSize = ::sqrt(coeffs.size()/columnIndices.size());

    // Create coefficient matrix: we must copy coeffs from CSR storage
    // Maybe it is possible to use CSR storage directly.
    
    // Colm- reserving 6*6 = 36 more spaces for coefficients
    std::vector< Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(36*(d.size() + l.size() + u.size() + 1));
    //-----------------------------------------------------------------------------
    //                  diagonal
    //-----------------------------------------------------------------------------
    label globalRowI = 0;
    forAll(d, cellI)
    {
        const scalarSquareMatrix& curD = d[cellI];

        for (label localRowI = 0; localRowI < 6; localRowI++)
        {
            for (label localColI = 0; localColI < 6; localColI++)
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
    // Colm: Add rigid-body block: 6x6 identity
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

//-----------------------------------------------------------------------------
//                  off-diagonal
//-----------------------------------------------------------------------------

    forAll(u, faceI)
    {
        const scalarSquareMatrix& curU = u[faceI];
        const scalarSquareMatrix& curL = l[faceI];

        const label globalRowI = 6*own[faceI];
        const label globalColI = 6*nei[faceI];

        //upper
        for (label localRowI = 0; localRowI < 6; localRowI++)
        {
            for (label localColI = 0; localColI < 6; localColI++)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalRowI + localRowI,
                        globalColI + localColI,
                        curU(localRowI, localColI)
                    )
                );
            }
        }
        //lower
        for (label localRowI = 0; localRowI < 6; localRowI++)
        {
            for (label localColI = 0; localColI < 6; localColI++)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalColI + localRowI,
                        globalRowI + localColI,
                        curL(localRowI, localColI)
                    )
                );
            }
        }


    }
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------

        // Insert triplets into the matrix
        //label bnRows = rowPointers.size()-1;
        //nRows = blockSize*bnRows;
    A.resize(nRows, nRows);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Compressing matrix is meant to help performance
    A.makeCompressed();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::scalar Foam::BlockEigenSolverOF::solve
(
    Foam::Field<Foam::scalarRectangularMatrix>& foamX,
    const Foam::Field<Foam::scalarRectangularMatrix>& foamB,
    Foam::scalarRectangularMatrix& sixDOFx
 )
{
    // Allow to run in parallel, where each core solves its own independent
    // problem: this will not be correct if the beam is split across cores!
    // if (Pstream::parRun())
    // {
    //     FatalErrorIn
    //     (
    //         "bool Foam::BlockEigenSolverOF::solve"
    //         "(...)"
    //     )   << "Eigen direct linear solver may not be run in parallel"
    //         << abort(FatalError);
    // }

    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A; // initialized in convertFoamMatrixToEigenMatrix funtion
    convertFoamMatrixToEigenMatrix(d_, l_, u_, own_, nei_, A);
    // Create Eigen source and solution vector from foam vectors
    //Eigen::Matrix<scalar, Eigen::Dynamic, 1> b;
    //Eigen::Matrix<scalar, Eigen::Dynamic, 1> x;

    const label nRows = A.rows();
    Info << "Number of rows in matrix called by solver = " << nRows << endl; 
    
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

    // Colm: rigidStart counter for starting work on rigidBody
    label rigidStart = nRows - 6;

    // Colm: defining variables for rigid body RHS
    // Colm: in future variables will be read in

    scalar rb_newmarkB_gamma = 0.5;
    scalar rb_newmarkB_beta = 0.25;
    scalar deltaT = 1;
    
    vector rb_v_old(1,1,1);
    vector rb_disp_old(1,1,1);
    vector rb_a(0.5,0.5,0.5);
    vector rb_a_old(0.5,0.5,0.5);

    vector rb_angMomentum_old(1,1,1);
    vector rb_tau(0.5,0.5,0.5);
    vector rb_tau_old(0.5,0.5,0.5);

    
    // Colm: Rigid body RHS

    // Colm: Linear Motion- displacement
    for (label i = 0; i < 3; ++i)
      {
	b(rigidStart + i) =
	  rb_disp_old[i]
	  + deltaT*rb_v_old[i]
	  + deltaT*deltaT*
	  (rb_newmarkB_beta*rb_a[i]
	   + (0.5 - rb_newmarkB_beta)*rb_a_old[i]);
      }

    // Colm: Angular Motion- correction to rotation
    for (label i = 0; i < 3; ++i)
      {
	b(rigidStart + 3 + i) =
	  deltaT*rb_angMomentum_old[i]
	  + deltaT*deltaT*
	  (rb_newmarkB_beta*rb_tau[i]
	  + (0.5 - rb_newmarkB_beta)*rb_tau_old[i]);
      }


    // Copy solution vector into Eigen vector
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
    Info << "Input sixDOFx:" << endl;
    Info << "______________" << endl;
    for (label i = 0; i < 6; ++i)
    {
      x(rigidStart + i) = 0.0;
      Info << "    " << x(rigidStart + i) << endl;
    }
    Info << "______________" << endl;


    // Calculate initial residual
    // Colm: +1 for extra 6 x 6 block
    const label nCells = d_.size() + 1;
    scalar initialResidual = 0;
    {

        Eigen::Matrix<scalar, Eigen::Dynamic, 1> Ax(nRows);
        Ax = A*x;

        Field<scalarRectangularMatrix> foamAx
        (
            nCells, scalarRectangularMatrix(6, 1, 0.0)
        );

        // Convert product
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<6; j++)
            {
                foamAx[i](j,0) = Ax[k++];
            }
        }

        Field<scalarRectangularMatrix> blockR(foamB - foamAx);

        //scalarRectangularMatrix norm(6, 1, 1.0); //this->normFactor(U, blockB);
        //initialResidual = cmptDivide(gSum(cmptMag(blockR)), norm);
        // initialResidual = sqrt(gSum(magSqr(blockR)));
        // Note: use sum instead of gSum
        initialResidual = sqrt(sum(magSqr(blockR)));
    }


    typedef enum {SparseLU, BiCGSTAB, GMRES, DGMRES, MINRES} solvers;
    solvers sol = SparseLU;

    switch(sol)
    {
        case SparseLU:
        {
            Eigen::SparseLU
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::COLAMDOrdering<int>
            > solver(A);
            x = solver.solve(b);
            //solverPerf.nIterations()++;
            break;
        }
        case BiCGSTAB:
        {
            Eigen::BiCGSTAB
            <
                Eigen::SparseMatrix<scalar>,
                // Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        case GMRES:
        {
            Eigen::GMRES
            <
                Eigen::SparseMatrix<scalar>,
                // Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        case DGMRES:
        {
            Eigen::DGMRES
            <
                Eigen::SparseMatrix<scalar>,
                // Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        case MINRES:
        {
            Eigen::MINRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::Lower|Eigen::Upper,
                Eigen::IdentityPreconditioner
                // Eigen::DiagonalPreconditioner<scalar>
                // Eigen::IncompleteLUT<scalar>
            > solver;
            solver.compute(A);
            x = solver.solve(b);
            //solverPerf.nIterations() = solver.iterations();
            break;
        }
        default:
            FatalErrorIn
            (
                "Foam::BlockSolverPerformance<Foam::vector6>"
                "Foam::BlockEigenSolverOF::solve"
                "("
                "    Field<Foam::vector6>& U,"
                "    const Field<Foam::vector6>& blockB"
                ")"
            )   << "Undefined Eugen solver."
                << abort(FatalError);
    }

    // We copy the results from the std::vector into the geometric field
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

    Info << "Output sixDOFx:" << endl;
    Info << "______________" << endl;
    for (label i = 0; i < 6; ++i)
      {
	sixDOFx(i,0) = x(rigidStart + i);
	Info << "    " << x(rigidStart + i) << endl;
      }
    Info << "______________" << endl;

    //

    // Calculate final residual
    // {
    //     Eigen::Matrix<scalar, Eigen::Dynamic, 1> p(nRows);
    //     p = A*x;

    //     Field<vector6> blockP(U.size());

    //     // Convert poroduct
    //     label k = 0;
    //     for (label i=0; i<nCells; i++)
    //     {
    //         for (label j=0; j<blockSize; j++)
    //         {
    //             blockP[i](j) = p[k++];
    //         }
    //     }

    //     Field<vector6> blockR(blockB - blockP);

    //     vector6 norm = vector6::one; //this->normFactor(U, blockB);

    //     solverPerf.finalResidual() =
    //         cmptDivide(gSum(cmptMag(blockR)), norm);
    // }


    return initialResidual;
}
