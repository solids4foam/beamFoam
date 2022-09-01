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

#include "BlockEigenSolver.H"

#include "OFstream.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "beamModel.H"

#include "multibeamFvBlockMatrix.H"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSolver, 0);

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockEigenSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockEigenSolver, asymMatrix
    );
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

void Foam::BlockEigenSolver::convertFoamMatrixToEigenMatrix
(
    const BlockLduMatrix<vector6>& matrix,
    Eigen::SparseMatrix<scalar>& A,
    Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b,
    Eigen::Matrix<scalar, Eigen::Dynamic, 1>& x
)
{
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName
            << ": copying matrix coefficients into Eigen format"
            << endl;
    }
    
    multibeamFvBlockMatrix& mbMatrix =
        const_cast<multibeamFvBlockMatrix&>
        (
            static_cast<const multibeamFvBlockMatrix&>(matrix_)
        );

    label nRows = 0;
        
    if (mbMatrix.hybrid())
    {
        // CSR matrix storage
        
        const labelList& rowPointers = mbMatrix.rowPointers();
        const labelList& columnIndices = mbMatrix.columnIndices();
        const scalarField& coeffs = mbMatrix.coeffs();

        // Create coefficient matrix: we must copy coeffs from CSR storage
        // Maybe it is possible to use CSR storage directly.
        std::vector< Eigen::Triplet<scalar> > coefficients;
        coefficients.reserve(coeffs.size());
    
        for (label rowI=1; rowI<rowPointers.size(); rowI++)
        {
            label colStart = rowPointers[rowI-1];
            label colEnd = rowPointers[rowI];

            for (label colI=colStart; colI<colEnd; colI++)
            {
                label i = rowI-1;
                label j = columnIndices[colI];
            
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>(i, j, coeffs[colI])
                );
            }
        }

        // Insert triplets into the matrix
        nRows = rowPointers.size()-1;
        A.resize(nRows, nRows);
        A.setFromTriplets(coefficients.begin(), coefficients.end());
    }
    else
    {
        // Block CSR matrix storage
        
        const labelList& rowPointers = mbMatrix.blockRowPointers();
        const labelList& columnIndices = mbMatrix.blockColumnIndices();
        const scalarField& coeffs = mbMatrix.blockCoeffs();

        label blockSize = ::sqrt(coeffs.size()/columnIndices.size());

        // Create coefficient matrix: we must copy coeffs from CSR storage
        // Maybe it is possible to use CSR storage directly.
        std::vector< Eigen::Triplet<scalar> > coefficients;
        coefficients.reserve(coeffs.size());

        label k=0;
        for (label rowI=1; rowI<rowPointers.size(); rowI++)
        {
            label colStart = rowPointers[rowI-1];
            label colEnd = rowPointers[rowI];

            for (label colI=colStart; colI<colEnd; colI++)
            {
                label i = blockSize*(rowI-1);
            
                for (label rI=0; rI<blockSize; rI++)
                { 
                    label j = blockSize*columnIndices[colI];
                    
                    for (label cI=0; cI<blockSize; cI++)
                    {
                        coefficients.push_back
                        (
                            Eigen::Triplet<scalar>(i, j++, coeffs[k++])
                        );
                    }
                    i++;
                }
            }
        }
        
        // Insert triplets into the matrix
        label bnRows = rowPointers.size()-1;
        nRows = blockSize*bnRows;
        A.resize(nRows, nRows);
        A.setFromTriplets(coefficients.begin(), coefficients.end());
    }

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
Foam::BlockEigenSolver::BlockEigenSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector6>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector6>(fieldName, matrix, dict)
{}

// ************************************************************************* //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector6>
Foam::BlockEigenSolver::solve
(
    Field<Foam::vector6>& U,
    const Field<Foam::vector6>& blockB
)
{
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector6>"
            "Foam::BlockEigenSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        )   << "Eigen direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Prepare solver performance
    BlockSolverPerformance<vector6> solverPerf
    (
        typeName,
        this->fieldName()
    );


    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A; //(nRows, nRows);
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b;
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x;

    convertFoamMatrixToEigenMatrix(matrix_, A, b, x);

    label nRows = A.rows();
    label nCells = blockB.size();
    label blockSize = 6;

    // Copy source vector into Eigen vector
    // Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(nRows);
    // label index = 0;
    // forAll(blockB, cellI)
    // {
    //     b(index++) = blockB[cellI](0);
    //     b(index++) = blockB[cellI](1);
    //     b(index++) = blockB[cellI](2);
    //     b(index++) = blockB[cellI](3);
    //     b(index++) = blockB[cellI](4);
    //     b(index++) = blockB[cellI](5);
    // }

    // Copy solution vector into Eigen vector
    // Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(nRows);
    // index = 0;
    // forAll(U, cellI)
    // {
    //     x(index++) = U[cellI](0);
    //     x(index++) = U[cellI](1);
    //     x(index++) = U[cellI](2);
    //     x(index++) = U[cellI](3);
    //     x(index++) = U[cellI](4);
    //     x(index++) = U[cellI](5);
    // }

    // Calculate initial residual
    {
        Eigen::Matrix<scalar, Eigen::Dynamic, 1> p(nRows);
        p = A*x;

        Field<vector6> blockP(nCells);
        
        // Convert poroduct
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<blockSize; j++)
            {
                blockP[i](j) = p[k++];
            }
        }

        Field<vector6> blockR(blockB - blockP);

        vector6 norm = vector6::one; //this->normFactor(U, blockB);

        solverPerf.initialResidual() =
            cmptDivide(gSum(cmptMag(blockR)), norm);
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
            solverPerf.nIterations()++;
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
            solverPerf.nIterations() = solver.iterations();
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
            solverPerf.nIterations() = solver.iterations();
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
            solverPerf.nIterations() = solver.iterations();
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
            solverPerf.nIterations() = solver.iterations();
            break;
        }
        default:
            FatalErrorIn
            (
                "Foam::BlockSolverPerformance<Foam::vector6>"
                "Foam::BlockEigenSolver::solve"
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
        U[cellI](0) = x(index++);
        U[cellI](1) = x(index++);
        U[cellI](2) = x(index++);
        U[cellI](3) = x(index++);
        U[cellI](4) = x(index++);
        U[cellI](5) = x(index++);
    }

    //
    if (nRows > nCells*blockSize)
    {
        multibeamFvBlockMatrix& mbMatrix =
            const_cast<multibeamFvBlockMatrix&>
            (
                static_cast<const multibeamFvBlockMatrix&>(matrix_)
            );

        for (label i=0; i<nRows; i++)
        {
            mbMatrix.solution()[i] = x(i);
        }
    }
    
    // Calculate final residual
    {
        Eigen::Matrix<scalar, Eigen::Dynamic, 1> p(nRows);
        p = A*x;

        Field<vector6> blockP(U.size());

        // Convert poroduct
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<blockSize; j++)
            {
                blockP[i](j) = p[k++];
            }
        }

        Field<vector6> blockR(blockB - blockP);

        vector6 norm = vector6::one; //this->normFactor(U, blockB);

        solverPerf.finalResidual() =
            cmptDivide(gSum(cmptMag(blockR)), norm);
    }
    
    return solverPerf;
}


