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

#include "blockLduSolvers.H"
#include "BlockEigenSparseLUSolver.H"
#include "addToRunTimeSelectionTable.H"
#include <Eigen/Sparse>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSparseLUSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockEigenSparseLUSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockEigenSparseLUSolver, asymMatrix
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenSparseLUSolver::BlockEigenSparseLUSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector6>& matrix,
    const dictionary& dict
)
:
    BlockEigenSolver(fieldName, matrix, dict)
{
    Info<< type() << " Eigen linear solver" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector6>
Foam::BlockEigenSparseLUSolver::solve
(
    Field<Foam::vector6>& U,
    const Field<Foam::vector6>& blockB
)
{
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector>"
            "Foam::BlockEigenSparseLUSolver::solve"
            "("
            "    Field<Foam::vector>& U,"
            "    const Field<Foam::vector>& blockB"
            ")"
        )   << "SparesLU direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();

    // Calculate te number of degrees of freedom
    const int m = calcDegreesOfFreedom(matrix_, twoD);

    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A(m, m);
    
    // Convert foam matrix to the Eigen matrix format
    //convertFoamMatrixToEigenMatrix(d, u, l, upperAddr, lowerAddr, twoD, A);
    convertFoamMatrixToEigenMatrix(matrix_, A);

    // std::cout << A;
    
    // Copy source vector into Eigen vector
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(m);
    label index = 0;
    forAll(blockB, rowI)
    {
        b(index++) = blockB[rowI](0);
        b(index++) = blockB[rowI](1);
        b(index++) = blockB[rowI](2);
        b(index++) = blockB[rowI](3);
        b(index++) = blockB[rowI](4);
        b(index++) = blockB[rowI](5);

        // if (!twoD)
        // {
        //     b(index++) = blockB[rowI].z();
        // }
    }

    // Copy source vector into Eigen vector
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(m);
    index = 0;
    forAll(U, rowI)
    {
        x(index++) = U[rowI](0);
        x(index++) = U[rowI](1);
        x(index++) = U[rowI](2);
        x(index++) = U[rowI](3);
        x(index++) = U[rowI](4);
        x(index++) = U[rowI](5);
    }

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> r = A*x-b;

    Field<Foam::vector6> R(U.size(), Foam::vector6::zero);
    index = 0;
    forAll(R, cellI)
    {
        R[cellI](0) = mag(r(index++));
        R[cellI](1) = mag(r(index++));
        R[cellI](2) = mag(r(index++));
        R[cellI](3) = mag(r(index++));
        R[cellI](4) = mag(r(index++));
        R[cellI](5) = mag(r(index++));
    }

    // Info << mag(R) << endl;
    
    // Optionally export system to Matlab
    if (writeMatlabFiles())
    {
        writeLinearSystemToMatlabFiles(A, b);
    }

    // Compressing matrix is meant to help performance
    A.makeCompressed();

    // Create scaling
    //Eigen::IterScaling<Eigen::SparseMatrix<double> >* scalPtr = NULL;

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": direct solution of sparse system" << endl;
    }

    // Now, solve the equilibrated linear system with any available solver
    // In this case we use a Direct sparse LU solver
    // Trial an error suggests COLAMD ordering is by far the best
    //Eigen::SimplicialLDLT
    Eigen::SparseLU
    <
        Eigen::SparseMatrix<scalar>, Eigen::COLAMDOrdering<int>
        // Eigen::SparseMatrix<scalar>, Eigen::AMDOrdering<int>
        // Eigen::SparseMatrix<scalar>, Eigen::MetisOrdering<int>
    > solver(A);

    // A.colPivHouseholderQr().solve(b);
    
    x = solver.solve(b);
    // Eigen::Matrix<scalar, Eigen::Dynamic, 1> x = solver.solve(b);

    // We copy the results from the std::vector into the geometric field
    // This can be costly
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": copying results into foam format" << nl
            << endl;
    }
    index = 0;
    forAll(U, cellI)
    {
        U[cellI](0) = x(index++);
        U[cellI](1) = x(index++);
        U[cellI](2) = x(index++);
        U[cellI](3) = x(index++);
        U[cellI](4) = x(index++);
        U[cellI](5) = x(index++);

        // if (!twoD)
        // {
        //     U[cellI].z() = x(index++);
        // }
    }

    /*
    // Testing: eigenproblem
    {
        Info<< nl << "Computing the eigenvalues and eigenvectors" << endl;

        // Only for symmetric matrices; only uses the lower
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(A);

        Info<< "Copying one of the eigenvalues and eigenvectors" << endl;
        index = 0;
        forAll(U, cellI)
        {
            U[cellI].x() = es.eigenvectors().col(m - 3)(index++);
            U[cellI].y() = es.eigenvectors().col(m - 3)(index++);

            if (!twoD)
            {
                U[cellI].z() =
                    es.eigenvectors().col(m - 3)(index++);
            }
        }

        // These methods are much slower and I am not sure if I trust their
        // answers
        //Eigen::JacobiSVD<Eigen::MatrixXf> svd
        // Eigen::BDCSVD<Eigen::MatrixXf> svd
        // (
        //     A, Eigen::ComputeThinU | Eigen::ComputeThinV
        // );

        // Info<< "Copying the 2nd last eigenvalues and eigenvectors" << endl;
        // index = 0;
        // forAll(U, cellI)
        // {
        //     U[cellI].x() = svd.matrixU()(index++, m - 3);
        //     U[cellI].y() = svd.matrixU()(index++, m - 3);

        //     if (!twoD)
        //     {
        //         U[cellI].z() =
        //             svd.matrixU()(index++, m - 3);
        //     }
        // }
    }*/


    return BlockSolverPerformance<Foam::vector6>
    (
        this->typeName,
        this->fieldName(),
        pTraits<Foam::vector6>::one,
        pTraits<Foam::vector6>::zero,
        0,
        true,
        false
    );
}


// ************************************************************************* //
