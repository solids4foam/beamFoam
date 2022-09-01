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

#include "BlockMklPardisoSolver.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "beamModel.H"

#include "multibeamFvBlockMatrix.H"

#include "mkl.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockMklPardisoSolver, 0);

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockMklPardisoSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockMklPardisoSolver, asymMatrix
    );
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

void Foam::BlockMklPardisoSolver::convertFoamMatrixToMklPardisoMatrix
(
    const BlockLduMatrix<vector6>& matrix,
    labelList& rowPointers,
    labelList& columnIndices,
    scalarField& coeffs
)
{
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName
            << ": copying matrix coefficients into MklPardiso format"
            << endl;
    }

    // Matrix addressing
    const lduAddressing& lduAddr = matrix.lduAddr();

    const unallocLabelList& losortAddr = lduAddr.losortAddr();
    const unallocLabelList& losortStartAddr = lduAddr.losortStartAddr();
    const unallocLabelList& ownerStartAddr = lduAddr.ownerStartAddr();
    // const unallocLabelList& upperAddr = lduAddr.upperAddr();
    // const unallocLabelList& lowerAddr = lduAddr.lowerAddr();

    // Grab matrix diagonal and off-diagonals
    const Field<tensor6>& d = matrix.diag().asSquare();
    const Field<tensor6>& l = matrix.lower().asSquare();
    const Field<tensor6>& u = matrix.upper().asSquare();

    const fvMesh& mesh
    (
        refCast<const fvMesh>
        (
            matrix.mesh()
        )
    );

    // const multibeamFvBlockMatrix& mbMatrix =
    //     static_cast<const multibeamFvBlockMatrix&>(matrix_);

    // Info << matrix_.type() << mbMatrix.type() << endl;
    
    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>
        (
            "beamProperties"
        );

    label blockSize = 6*6;

    rowPointers = beam.csrAddr().rowPointers();
    columnIndices = beam.csrAddr().columnIndices();

    coeffs.setSize(blockSize*columnIndices.size(), 0);
    label coeffI = 0;

    for (label rowI=1; rowI<rowPointers.size(); rowI++)
    {
        // Lower
        for
        (
            label j=losortStartAddr[rowI-1];
            j<losortStartAddr[rowI];
            j++
        )
        {
            // label nei = rowI-1;
            // label own = losortAddr[j];

            // label faceI = -1;

            // for
            // (
            //     label fI=ownerStartAddr[own];
            //     fI<ownerStartAddr[own+1];
            //     fI++
            // )
            // {
            //     if (upperAddr[fI] == nei)
            //     {
            //         faceI = fI;
            //         break;
            //     }
            // }

            label faceI = losortAddr[j];
                
            const tensor6& curCoeff = l[faceI];

            // for (label jj=0; jj<6; jj++)
            // {
            //     for (label ii=0; ii<6; ii++)
            //     {
            //         // column-major order for 1-based indexing
            //         coeffs[coeffI++] = curCoeff(ii,jj);
            //     }
            // }
            for (label ii=0; ii<6; ii++)
            {
                for (label jj=0; jj<6; jj++)
                {
                    // row-major order for 0-based indexing
                    coeffs[coeffI] = curCoeff(ii,jj);
                    coeffI++;
                }
            }
        }

        // Diagonal
        {
            const tensor6& curCoeff = d[rowI-1];

            // for (label jj=0; jj<6; jj++)
            // {
            //     for (label ii=0; ii<6; ii++)
            //     {
            //         // column-major order for 1-based indexing
            //         coeffs[coeffI++] = curCoeff(ii,jj);
            //     }
            // }
            for (label ii=0; ii<6; ii++)
            {
                for (label jj=0; jj<6; jj++)
                {
                    // row-major order for 0-based indexing
                    coeffs[coeffI] = curCoeff(ii,jj);
                    coeffI++;
                }
            }
        }

        // Upper
        for
        (
            label j=ownerStartAddr[rowI-1];
            j<ownerStartAddr[rowI];
            j++
        )
        {
            label faceI = j;
            
            const tensor6& curCoeff = u[faceI];

            // for (label jj=0; jj<6; jj++)
            // {
            //     for (label ii=0; ii<6; ii++)
            //     {
            //         // column-major order for 1-based indexing
            //         coeffs[coeffI++] = curCoeff(ii,jj);
            //     }
            // }
            for (label ii=0; ii<6; ii++)
            {
                for (label jj=0; jj<6; jj++)
                {
                    // row-major order for 0-based indexing
                    coeffs[coeffI] = curCoeff(ii,jj);
                    coeffI++;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockMklPardisoSolver::BlockMklPardisoSolver
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
Foam::BlockMklPardisoSolver::solve
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
            "Foam::BlockMklPardisoSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        )   << "MklPardiso direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Prepare solver performance
    BlockSolverPerformance<vector6> solverPerf
    (
        typeName,
        this->fieldName()
    );
    
    vector6 norm = vector6::one; //this->normFactor(U, blockB);

    // Calculate initial residual
    {
        Field<vector6> p(U.size());
        matrix_.Amul(p, U);
        Field<vector6> r(blockB - p);

        solverPerf.initialResidual() = cmptDivide(gSum(cmptMag(r)), norm);
    }

    labelList rowPointers;
    labelList columnIndices;
    scalarField coeffs;

    convertFoamMatrixToMklPardisoMatrix
    (
        matrix_,
        rowPointers,
        columnIndices,
        coeffs
    );

    // // Convert to one based indexing
    // forAll(rowPointers, i)
    // {
    //     rowPointers[i] += 1;
    // }
    // forAll(columnIndices, i)
    // {
    //     columnIndices[i] += 1;
    // }
    
    // Block size size
    MKL_INT blockSize = 6;

    // Number of cells
    MKL_INT nCells = U.size();

    // Matrix size (block matrix)
    MKL_INT n = nCells;

    // RHS and solution vectors contruction
    scalarField b(blockSize*n, 0);
    scalarField x(blockSize*n, 0);

    // RHS and solution vector conversion
    label k = 0;
    for (label i=0; i<nCells; i++)
    {
        for (label j=0; j<blockSize; j++)
        {
            b[k] = blockB[i](j);
            x[k] = U[i](j);
            k++;
        }
    }

    // // Calculate initial residual
    // {
    //     scalarField p(x.size());
    //     MatVecMul(rowPointers, columnIndices, coeffs, x, p);
    //     scalarField r(b - p);

    //     scalar res = ::sqrt(sum(magSqr(r)))/::sqrt(sum(magSqr(b)));
    //     solverPerf.initialResidual() = vector6::zero;
    //     solverPerf.initialResidual()(0) = res;
        
    //     // solverPerf.initialResidual() = cmptDivide(gSum(cmptMag(rb)), norm);
    //     // solverPerf.nIterations()++;
    // }
    
    // Matrix type: real unsymmetric
    MKL_INT mtype = 11;

    // Number of right hand sides.
    MKL_INT nrhs = 1;

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void* pt[64];

    // Pardiso control parameters
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;

    // Auxiliary variables
    scalar ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */

    // Set number of threads
    mkl_set_num_threads(1);

    // Setup Pardiso control parameters
    for ( label i=0; i<64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;          /* No solver default */
    iparm[1] = 3;          /* Fill-in reordering from METIS */
    iparm[3] = 0;          /* No iterative-direct algorithm */
    iparm[4] = 0;          /* No user fill-in reducing permutation */
    iparm[5] = 0;          /* Write solution into x */
    iparm[6] = 0;          /* Not in use */
    iparm[7] = 1;          /* Max numbers of iterative refinement steps 2 zt*/
    iparm[8] = 0;          /* Not in use */
    iparm[9] = 13;         /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 0;         /* Use nonsymmetric permutation and scaling MPS 1 zt*/
    iparm[11] = 0;         /* Conjugate transposed/transpose solve */
    iparm[12] = 0;         /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;         /* Output: Number of perturbed pivots */
    iparm[14] = 0;         /* Not in use */
    iparm[15] = 0;         /* Not in use */
    iparm[16] = 0;         /* Not in use */
    iparm[17] = -1;        /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;        /* Output: Mflops for LU factorization */
    iparm[19] = 0;         /* Output: Numbers of CG Iterations */
    iparm[23] = 0;         /* Parallel factorization control */
    iparm[26] = 0;         /* Matrix checker */
    iparm[36] = blockSize; /* Block size */
    iparm[34] = 1;          /* Zero-based indexing of columns and rows */
    iparm[59] = 0;         /* IC/OOC */
    maxfct = 1;            /* Maximum number of numerical factorizations. */
    mnum = 1;              /* Which factorization to use. */
    msglvl = 0;            /* Print statistical information  */
    error = 0;             /* Initialize error flag */

    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver
    for (label i=0; i<64; i++)
    {
        pt[i] = 0;
    }

    // Reordering and Symbolic Factorization. This step also allocates
    // all memory that is necessary for the factorization
    phase = 11;
    PARDISO
    (
        pt,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        coeffs.cdata(),
        rowPointers.cdata(),
        columnIndices.cdata(),
        &idum,
        &nrhs,
        iparm,
        &msglvl,
        &ddum,
        &ddum,
        &error
    );

    if (error!=0)
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector6>"
            "Foam::BlockMklPardisoSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        ) << "ERROR during symbolic factorization: " << error
          << abort(FatalError);
    }
    // Info << "Reordering completed ... " << endl;
    // Info << "Number of nonzeros in factors = " << iparm[17] << endl;
    // Info << "Number of factorization MFLOPS = " << iparm[18] << endl;

    // Numerical factorization
    phase = 22;
    PARDISO
    (
        pt,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        coeffs.cdata(),
        rowPointers.cdata(),
        columnIndices.cdata(),
        &idum,
        &nrhs,
        iparm,
        &msglvl,
        &ddum,
        &ddum,
        &error
    );
    
    if (error != 0)
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector6>"
            "Foam::BlockMklPardisoSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        ) << "ERROR during numerical factorization: " << error
          << abort(FatalError);
    }
    // Info << "Factorization completed ... " << endl;

    // Solve
    phase = 33;
    PARDISO
    (
        pt,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        coeffs.cdata(),
        rowPointers.cdata(),
        columnIndices.cdata(),
        &idum,
        &nrhs,
        iparm,
        &msglvl,
        b.data(),
        x.data(),
        &error
    );

    if (error!=0)
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector6>"
            "Foam::BlockMklPardisoSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        ) << "ERROR during solution: " << error
          << abort(FatalError);
    }

    // Convert solution vector
    k = 0;
    for (label i=0; i<nCells; i++)
    {
        for (label j=0; j<blockSize; j++)
        {
            U[i](j) = x[k];
            k++;
        }
    }

    // Calculate final residual
    {
        Field<vector6> p(U.size());
        matrix_.Amul(p, U);
        Field<vector6> r(blockB - p);

        solverPerf.finalResidual() = cmptDivide(gSum(cmptMag(r)), norm);
        solverPerf.nIterations()++;
    }

    // // Calculate final residual
    // {
    //     scalarField p(x.size());
    //     MatVecMul(rowPointers, columnIndices, coeffs, x, p);
    //     scalarField r(b - p);

    //     scalar res = sqrt(sum(magSqr(r)))/sqrt(sum(magSqr(b)));
    //     solverPerf.finalResidual() = vector6::zero;
    //     solverPerf.finalResidual()(0) = res;
        
    //     // solverPerf.finalResidual() = cmptDivide(gSum(cmptMag(rb)), norm);
    //     solverPerf.nIterations()++;
    // }
    
    // Termination and release of memory
    phase = -1; // Release internal memory
    PARDISO
    (
        pt,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        &ddum,
        &idum, //rowPointers.cdata(),
        &idum, //columnIndices.cdata(),
        &idum,
        &nrhs,
        iparm,
        &msglvl,
        &ddum,
        &ddum,
        &error
    );

    return solverPerf;
}


void Foam::BlockMklPardisoSolver::MatVecMul
(
    const labelList& rowPointers,
    const labelList& columnIndices,
    const scalarField& coeffs,
    const scalarField& x,
    scalarField& p
) const
{
    label blockSize = ::sqrt(coeffs.size()/columnIndices.size());
    // Info << "blockSize = " << blockSize << endl;

    label nRows = rowPointers.size()-1;

    label k=0;
    for (label i=0; i<nRows; i++)
    {
        const label jStart = rowPointers[i];
        const label jEnd = rowPointers[i+1];

        for (label j=jStart; j<jEnd; j++)
        {
            const label cI = columnIndices[j];

            for (label ib=blockSize*i; ib<blockSize*(i+1); ib++)
            {
                for (label jb=blockSize*cI; jb<blockSize*(cI+1); jb++)
                {
                    p[ib] += coeffs[k]*x[jb];
                    k++;
                }
            }
        }
    }
}
