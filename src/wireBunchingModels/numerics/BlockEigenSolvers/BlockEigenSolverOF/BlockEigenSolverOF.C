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
    const Field<scalarSquareMatrix>& d_,
    const Field<scalarSquareMatrix>& l_,
    const Field<scalarSquareMatrix>& u_,
    const labelList& own,
    const labelList& nei,
    Eigen::SparseMatrix<scalar>& A
)
{
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName
            << ": copying matrix coefficients into Eigen format"
            << endl;
    }

    //label nRows = 6*d.size();

    // Block CSR matrix storage

    //const labelList& rowPointers = mbMatrix.blockRowPointers();
    //const labelList& columnIndices = mbMatrix.blockColumnIndices();
    //const scalarField& coeffs = mbMatrix.blockCoeffs();

    //label blockSize = ::sqrt(coeffs.size()/columnIndices.size());

    // Create coefficient matrix: we must copy coeffs from CSR storage
    // Maybe it is possible to use CSR storage directly.
    std::vector< Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(coeffs.size());
    //-----------------------------------------------------------------------------
    //                  diagonal
    //-----------------------------------------------------------------------------
    label globalRowI = 0;
    forAll(d, cellI)
    {
        scalarSquareMatrix& curD = d[cellI];

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
//-----------------------------------------------------------------------------
//                  off-diagonal
//-----------------------------------------------------------------------------
    globalRowI = 0;
    forAll(u, faceI)
    {
        scalarSquareMatrix& curU = u[faceI];
        scalarSquareMatrix& curL = l[faceI];
        const label owner = own[faceI];
        const label neighbour = nei[faceI];


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
    //-----------------------------------------------------------------------------
    //-----------------------------------------------------------------------------

        // Insert triplets into the matrix
        label bnRows = rowPointers.size()-1;
        nRows = blockSize*bnRows;
        A.resize(nRows, nRows);
        A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Compressing matrix is meant to help performance
    A.makeCompressed();

    // Set rhs and solution vectors
    b.resize(nRows);
    x.resize(nRows);
    for (label i=0; i<nRows; i++)
    {
        b(i) = mbMatrix.rhs()[i];
        x(i) = mbMatrix.solution()[i];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenSolverOF::BlockEigenSolverOF
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei,
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
    const Foam::Field<Foam::scalarRectangularMatrix>& foamB
)
{
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "bool Foam::BlockEigenSolverOF::solve"
            "(...)"
        )   << "Eigen direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A; //(nRows, nRows);
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x;


    convertFoamMatrixToEigenMatrix(d, l, u, A);


    label nRows = A.rows();
    //label nCells = blockB.size();
    label nCells= foamB.size();
    label blockSize = 6;

    //Copy source vector into Eigen vector

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

    // Calculate initial residual
    scalar initialResidual = 0;
    {

        Eigen::Matrix<scalar, Eigen::Dynamic, 1> p(nRows);
        p = A*x;

        //Field<vector6> blockP(nCells);
        Field<scalarRectangularMatrix> blockP
    (
        nCells, scalarRectangularMatrix(6, 1, 0.0)
    );
        // Convert poroduct
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<blockSize; j++)
            {
                blockP[i](j,0) = p[k++];
            }
        }

        Field<scalarRectangularMatrix> blockR(blockB - blockP);

        scalarRectangularMatrix norm(6, 1, 1.0); //this->normFactor(U, blockB);

        initialResidual = cmptDivide(gSum(cmptMag(blockR)), norm);
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
    label index = 0;
    forAll(U, cellI)
    {
        foamX[cellI](0,0) = x(index++);
        foamX[cellI](1,0) = x(index++);
        foamX[cellI](2,0) = x(index++);
        foamX[cellI](3,0) = x(index++);
        foamX[cellI](4,0) = x(index++);
        foamX[cellI](5,0) = x(index++);
    }

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


