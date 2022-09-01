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

#include "BlockPetscSolver.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// #include "petscksp.h"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockPetscSolver, 0);

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockPetscSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockPetscSolver, asymMatrix
    );
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

bool Foam::BlockPetscSolver::checkTwoD() const
{
    return false;
}


int Foam::BlockPetscSolver::calcDegreesOfFreedom
(
    const BlockLduMatrix<vector6>& matrix,
    const bool twoD
) const
{
    // Number of vectors on the diagonal
    const int diagSize = matrix.diag().asSquare().size();

    if (twoD)
    {
        // Two components in each vector
        return 4*diagSize;
    }

    // Info << "diagSize: " << diagSize << endl;
    // Three components in each vector
    return 6*diagSize;
}


void Foam::BlockPetscSolver::convertFoamMatrixToPetscMatrix
(
    const BlockLduMatrix<vector6>& matrix,
    Mat& A
)
{
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName
            << ": copying matrix coefficients into Petsc format"
            << endl;
    }

    // Matrix addressing
    const unallocLabelList& lowerAddr = matrix.mesh().lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix.mesh().lduAddr().upperAddr();

    // Grab matrix diagonal and off-diagonals
    const Field<tensor6>& d = matrix.diag().asSquare();
    const Field<tensor6>& l = matrix.lower().asSquare();
    const Field<tensor6>& u = matrix.upper().asSquare();

    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();

    // Calculate the number of degrees of freedom
    const PetscInt n = calcDegreesOfFreedom(matrix, twoD);

    // Calculate the number of non-zero per row
    const PetscInt nz = 3*6;

    PetscErrorCode ierr;
    
    ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, n, n, nz, NULL, &A);
    CHKERRV(ierr);

    // Insert diagonal
    label ID = -1;
    forAll(d, dI)
    {
        if (twoD)
        {
            // 2-D
            ID = 4*dI;
        }
        else
        {
            // 3-D
            ID = 6*dI;
        }

        ierr = MatSetValue(A, ID, ID, d[dI](0,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID, ID+1, d[dI](0,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID, ID+2, d[dI](0,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID, ID+3, d[dI](0,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID, ID+4, d[dI](0,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID, ID+5, d[dI](0,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, ID+1, ID, d[dI](1,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+1, ID+1, d[dI](1,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+1, ID+2, d[dI](1,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+1, ID+3, d[dI](1,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+1, ID+4, d[dI](1,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+1, ID+5, d[dI](1,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, ID+2, ID, d[dI](2,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+2, ID+1, d[dI](2,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+2, ID+2, d[dI](2,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+2, ID+3, d[dI](2,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+2, ID+4, d[dI](2,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+2, ID+5, d[dI](2,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, ID+3, ID, d[dI](3,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+3, ID+1, d[dI](3,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+3, ID+2, d[dI](3,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+3, ID+3, d[dI](3,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+3, ID+4, d[dI](3,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+3, ID+5, d[dI](3,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, ID+4, ID, d[dI](4,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+4, ID+1, d[dI](4,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+4, ID+2, d[dI](4,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+4, ID+3, d[dI](4,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+4, ID+4, d[dI](4,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+4, ID+5, d[dI](4,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, ID+5, ID, d[dI](5,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+5, ID+1, d[dI](5,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+5, ID+2, d[dI](5,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+5, ID+3, d[dI](5,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+5, ID+4, d[dI](5,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, ID+5, ID+5, d[dI](5,5), INSERT_VALUES);
        CHKERRV(ierr);
    }
    
    // Insert off-diagonal
    label rowI = -1;
    label columnI = -1;
    forAll(u, uI)
    {
        const label own = lowerAddr[uI];
        const label nei = upperAddr[uI];

        const tensor6& upper = u[uI];
        const tensor6& lower = l[uI];

        if (twoD)
        {
            rowI = 4*own;
            columnI = 4*nei;
        }
        else
        {
            rowI = 6*own;
            columnI = 6*nei;
        }

        // Upper

        ierr = MatSetValue(A, rowI, columnI, upper(0,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI, columnI+1, upper(0,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI, columnI+2, upper(0,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI, columnI+3, upper(0,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI, columnI+4, upper(0,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI, columnI+5, upper(0,5), INSERT_VALUES);
        CHKERRV(ierr);
        
        ierr = MatSetValue(A, rowI+1, columnI, upper(1,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+1, columnI+1, upper(1,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+1, columnI+2, upper(1,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+1, columnI+3, upper(1,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+1, columnI+4, upper(1,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+1, columnI+5, upper(1,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, rowI+2, columnI, upper(2,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+2, columnI+1, upper(2,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+2, columnI+2, upper(2,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+2, columnI+3, upper(2,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+2, columnI+4, upper(2,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+2, columnI+5, upper(2,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, rowI+3, columnI, upper(3,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+3, columnI+1, upper(3,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+3, columnI+2, upper(3,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+3, columnI+3, upper(3,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+3, columnI+4, upper(3,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+3, columnI+5, upper(3,5), INSERT_VALUES);
        CHKERRV(ierr);
        
        ierr = MatSetValue(A, rowI+4, columnI, upper(4,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+4, columnI+1, upper(4,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+4, columnI+2, upper(4,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+4, columnI+3, upper(4,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+4, columnI+4, upper(4,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+4, columnI+5, upper(4,5), INSERT_VALUES);
        CHKERRV(ierr);
        
        ierr = MatSetValue(A, rowI+5, columnI, upper(5,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+5, columnI+1, upper(5,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+5, columnI+2, upper(5,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+5, columnI+3, upper(5,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+5, columnI+4, upper(5,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, rowI+5, columnI+5, upper(5,5), INSERT_VALUES);
        CHKERRV(ierr);

        // Lower
        
        ierr = MatSetValue(A, columnI, rowI, lower(0,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI, rowI+1, lower(0,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI, rowI+2, lower(0,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI, rowI+3, lower(0,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI, rowI+4, lower(0,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI, rowI+5, lower(0,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, columnI+1, rowI, lower(1,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+1, rowI+1, lower(1,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+1, rowI+2, lower(1,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+1, rowI+3, lower(1,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+1, rowI+4, lower(1,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+1, rowI+5, lower(1,5), INSERT_VALUES);
        CHKERRV(ierr);
        
        ierr = MatSetValue(A, columnI+2, rowI, lower(2,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+2, rowI+1, lower(2,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+2, rowI+2, lower(2,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+2, rowI+3, lower(2,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+2, rowI+4, lower(2,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+2, rowI+5, lower(2,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, columnI+3, rowI, lower(3,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+3, rowI+1, lower(3,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+3, rowI+2, lower(3,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+3, rowI+3, lower(3,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+3, rowI+4, lower(3,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+3, rowI+5, lower(3,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, columnI+4, rowI, lower(4,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+4, rowI+1, lower(4,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+4, rowI+2, lower(4,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+4, rowI+3, lower(4,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+4, rowI+4, lower(4,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+4, rowI+5, lower(4,5), INSERT_VALUES);
        CHKERRV(ierr);

        ierr = MatSetValue(A, columnI+5, rowI, lower(5,0), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+5, rowI+1, lower(5,1), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+5, rowI+2, lower(5,2), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+5, rowI+3, lower(5,3), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+5, rowI+4, lower(5,4), INSERT_VALUES);
        CHKERRV(ierr);
        ierr = MatSetValue(A, columnI+5, rowI+5, lower(5,5), INSERT_VALUES);
        CHKERRV(ierr);
    }

    
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockPetscSolver::BlockPetscSolver
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
Foam::BlockPetscSolver::solve
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
            "Foam::BlockPetscSolver::solve"
            "("
            "    Field<Foam::vector>& U,"
            "    const Field<Foam::vector>& blockB"
            ")"
        )   << "Petsc direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // PetscErrorCode ierr = PetscInitializeNoArguments();
    
    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();

    // Calculate the number of degrees of freedom
    const PetscInt n = calcDegreesOfFreedom(this->matrix_, twoD);
    
    // Create Eigen sparse matrix and set coeffs
    Mat A;
    
    // Convert foam matrix to the Petsc matrix format
    convertFoamMatrixToPetscMatrix(matrix_, A);

    // Copy source vector into Eigen vector
    Vec x, b, r;

    PetscErrorCode ierr;
    
    ierr = VecCreate(PETSC_COMM_WORLD, &x);
    CHKERRCONTINUE(ierr);
    ierr = VecSetType(x, VECSEQ);
    CHKERRCONTINUE(ierr);
    ierr = PetscObjectSetName((PetscObject) x, "Solution");
    CHKERRCONTINUE(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, n);
    CHKERRCONTINUE(ierr);
    
    ierr = VecDuplicate(x, &b);
    CHKERRCONTINUE(ierr);
    ierr = VecDuplicate(x, &r);
    CHKERRCONTINUE(ierr);
    
    // Copy source vector into Petsc vector
    label index = 0;
    forAll(blockB, rowI)
    {
        VecSetValue(b, index++, blockB[rowI](0), INSERT_VALUES);
        VecSetValue(b, index++, blockB[rowI](1), INSERT_VALUES);
        VecSetValue(b, index++, blockB[rowI](2), INSERT_VALUES);
        VecSetValue(b, index++, blockB[rowI](3), INSERT_VALUES);
        VecSetValue(b, index++, blockB[rowI](4), INSERT_VALUES);
        VecSetValue(b, index++, blockB[rowI](5), INSERT_VALUES);
    }

    // Copy solution vector into Petsc vector
    index = 0;
    forAll(U, rowI)
    {
        VecSetValue(x, index++, U[rowI](0), INSERT_VALUES);
        VecSetValue(x, index++, U[rowI](1), INSERT_VALUES);
        VecSetValue(x, index++, U[rowI](2), INSERT_VALUES);
        VecSetValue(x, index++, U[rowI](3), INSERT_VALUES);
        VecSetValue(x, index++, U[rowI](4), INSERT_VALUES);
        VecSetValue(x, index++, U[rowI](5), INSERT_VALUES);
    }

    // Assembling matrix
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRCONTINUE(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRCONTINUE(ierr);

    // Create linear solver
    KSP ksp;        // linear solver context
    PC pc;          // preconditioner context
    PetscReal norm; // norm of solution error
    PetscInt its;

    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRCONTINUE(ierr);

    // Set operators
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRCONTINUE(ierr);

    ierr = KSPGetPC(ksp, &pc);
    CHKERRCONTINUE(ierr);

    ierr = PCSetType(pc, PCILU); //PCHYPRE); //PCILU); //PCCHOLESKY);//PCILU);
    CHKERRCONTINUE(ierr);
    
    ierr = KSPSetTolerances
        (ksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    CHKERRCONTINUE(ierr);

    // Set KSP options
    ierr = KSPSetFromOptions(ksp);
    CHKERRCONTINUE(ierr);

    KSPSetType(ksp, KSPBCGS);
    
    // Solve the linear system
    ierr = KSPSolve(ksp, b, x);
    CHKERRCONTINUE(ierr);
    
    // View solver info
    ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
    CHKERRCONTINUE(ierr);

    // Check the solution and clean up
    ierr = MatMult(A, x, r);
    CHKERRCONTINUE(ierr);

    ierr = VecAXPY(r, -1.0, b);
    CHKERRCONTINUE(ierr);
    ierr = VecNorm(r, NORM_2, &norm);
    CHKERRCONTINUE(ierr);
    ierr = KSPGetIterationNumber(ksp, &its);
    CHKERRCONTINUE(ierr);
    ierr = PetscPrintf
    (
        PETSC_COMM_WORLD,
        "Norm of residual %g, Iterations %D\n",
        (double)norm,
        its
    );
    CHKERRCONTINUE(ierr);
    
    
    // We copy the results from the Petsc vector into the geometric field
    // This can be costly
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": copying results into foam format" << nl
            << endl;
    }
    index = 0;
    PetscScalar value;
    forAll(U, cellI)
    {
        ierr = VecGetValues(x, 1, &index, &value);
        CHKERRCONTINUE(ierr);
        U[cellI](0) = value;
        index++;
        
        ierr = VecGetValues(x, 1, &index, &value);
        CHKERRCONTINUE(ierr);
        U[cellI](1) = value;
        index++;
        
        ierr = VecGetValues(x, 1, &index, &value);
        CHKERRCONTINUE(ierr);
        U[cellI](2) = value;
        index++;
        
        ierr = VecGetValues(x, 1, &index, &value);
        CHKERRCONTINUE(ierr);
        U[cellI](3) = value;
        index++;
        
        ierr = VecGetValues(x, 1, &index, &value);
        CHKERRCONTINUE(ierr);
        U[cellI](4) = value;
        index++;
        
        ierr = VecGetValues(x, 1, &index, &value);
        CHKERRCONTINUE(ierr);
        U[cellI](5) = value;
        index++;
    }
    
    // // Eigen::Matrix<scalar, Eigen::Dynamic, 1> r = A*x-b;
    // double ir = (A*x - b).norm() / b.norm();
    
    // if (BlockLduSolver::debug)
    // {
    //     Info<< this->typeName << ": direct solution of sparse system" << endl;
    // }

    // // Now, solve the equilibrated linear system with any available solver
    // // In this case we use a Direct sparse LU solver
    // // Trial an error suggests COLAMD ordering is by far the best
    // //Eigen::SimplicialLDLT
    // Eigen::SparseLU
    // <
    //     Eigen::SparseMatrix<scalar>, Eigen::COLAMDOrdering<int>
    //     // Eigen::SparseMatrix<scalar>, Eigen::AMDOrdering<int>
    //     // Eigen::SparseMatrix<scalar>, Eigen::MetisOrdering<int>
    // > solver(A);

    // x = solver.solve(b);

    // double fr = (A*x - b).norm() / b.norm();
    // Info << "b.norm() = " << b.norm() << endl;
    
    // // We copy the results from the std::vector into the geometric field
    // // This can be costly
    // if (BlockLduSolver::debug)
    // {
    //     Info<< this->typeName << ": copying results into foam format" << nl
    //         << endl;
    // }
    // index = 0;
    // forAll(U, cellI)
    // {
    //     U[cellI](0) = x(index++);
    //     U[cellI](1) = x(index++);
    //     U[cellI](2) = x(index++);
    //     U[cellI](3) = x(index++);
    //     U[cellI](4) = x(index++);
    //     U[cellI](5) = x(index++);
    // }

    Foam::vector6 initRes = Foam::vector6::zero;
    initRes(0) = 1;
    // initRes(1) = ir;
    // initRes(2) = ir;
    // initRes(3) = ir;
    // initRes(4) = ir;
    // initRes(5) = ir;
    
    Foam::vector6 finalRes = Foam::vector6::zero;
    finalRes(0) = 0;
    // finalRes(1) = fr;
    // finalRes(2) = fr;
    // finalRes(3) = fr;
    // finalRes(4) = fr;
    // finalRes(5) = fr;
    
    return BlockSolverPerformance<Foam::vector6>
    (
        this->typeName,
        this->fieldName(),
        initRes,
        finalRes,
        its,
        true,
        false
    );
}
