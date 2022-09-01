/*---------------------------------------------------------------------------* \
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

#include "BlockMklClusterPardisoSolver.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "multibeamFvBlockMatrix.H"
#include "beamModel.H"


#include "mpi.h"
#include "mkl.h"
#include "mkl_cluster_sparse_solver.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockMklClusterPardisoSolver, 0);

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockMklClusterPardisoSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVector6Solver, BlockMklClusterPardisoSolver, asymMatrix
    );
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

void Foam::BlockMklClusterPardisoSolver::
convertFoamMatrixToMklClusterPardisoMatrix
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
            << ": copying matrix coefficients into MklClusterPardiso format"
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

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>
        (
            "beamProperties"
        );

    label blockSize = 6*6;

    rowPointers = beam.csrAddr().rowPointers();
    columnIndices = beam.csrAddr().columnIndices();

    label nCells = rowPointers.size()-1;

    coeffs.setSize(blockSize*columnIndices.size(), 0);
    label coeffI = 0;
    for (label rowI=1; rowI<rowPointers.size(); rowI++)
    {
        // Add lower nei processor cells
        if (beam.csrAddr().procCells().size())
        {
            if
            (
                beam.csrAddr().procCells()[rowI-1][0] != -1
             && beam.csrAddr().procCells()[rowI-1][2]
              < beam.csrAddr().globalNCellsOffset()
            )
            {
                label patchI = beam.csrAddr().procCells()[rowI-1][0];
                label faceI = beam.csrAddr().procCells()[rowI-1][1];
                
                const tensor6& curCoeff = 
                    this->matrix_.coupleUpper()[patchI].asSquare()[faceI];

                for (label ii=0; ii<6; ii++)
                {
                    for (label jj=0; jj<6; jj++)
                    {
                        // row-major order for 0-based indexing
                        coeffs[coeffI++] = -curCoeff(ii,jj);
                    }
                }
            }
        }
        
        // Lower
        for
        (
            label j=losortStartAddr[rowI-1];
            j<losortStartAddr[rowI];
            j++
        )
        {
            label faceI = losortAddr[j];
            
            const tensor6& curCoeff = l[faceI];

            for (label ii=0; ii<6; ii++)
            {
                for (label jj=0; jj<6; jj++)
                {
                    // row-major order for 0-based indexing
                    coeffs[coeffI++] = curCoeff(ii,jj);
                }
            }
        }

        // Diagonal
        {
            const tensor6& curCoeff = d[rowI-1];

            for (label ii=0; ii<6; ii++)
            {
                for (label jj=0; jj<6; jj++)
                {
                    // row-major order for 0-based indexing
                    coeffs[coeffI++] = curCoeff(ii,jj);
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

            for (label ii=0; ii<6; ii++)
            {
                for (label jj=0; jj<6; jj++)
                {
                    // row-major order for 0-based indexing
                    coeffs[coeffI++] = curCoeff(ii,jj);
                }
            }
        }

        // Add upper nei processor cells
        if (beam.csrAddr().procCells().size())
        {
            if
            (
                beam.csrAddr().procCells()[rowI-1][0] != -1
             && beam.csrAddr().procCells()[rowI-1][2]
             >= (beam.csrAddr().globalNCellsOffset() + nCells)
            )
            {
                label patchI = beam.csrAddr().procCells()[rowI-1][0];
                label faceI = beam.csrAddr().procCells()[rowI-1][1];

                const tensor6& curCoeff = 
                    this->matrix_.coupleUpper()[patchI].asSquare()[faceI];

                for (label ii=0; ii<6; ii++)
                {
                    for (label jj=0; jj<6; jj++)
                    {
                        // row-major order for 0-based indexing
                        coeffs[coeffI++] = -curCoeff(ii,jj);
                    }
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockMklClusterPardisoSolver::BlockMklClusterPardisoSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector6>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector6>(fieldName, matrix, dict),
    iparm40_(-1),
    iparm41_(-1),
    globalNCells_(0)
{}

// ************************************************************************* //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector6>
Foam::BlockMklClusterPardisoSolver::solve
(
    Field<Foam::vector6>& U,
    const Field<Foam::vector6>& blockB
)
{
    if (!Pstream::parRun())
    {
        return sequentialSolve(U, blockB);
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

    // labelList rowPointers;
    // labelList columnIndices;
    // scalarField coeffs;

    // convertFoamMatrixToMklClusterPardisoMatrix
    // (
    //     matrix_,
    //     rowPointers,
    //     columnIndices,
    //     coeffs
    // );

    multibeamFvBlockMatrix& mbMatrix =
        const_cast<multibeamFvBlockMatrix&>
        (
            static_cast<const multibeamFvBlockMatrix&>(matrix_)
        );

    const labelList& rowPointers = mbMatrix.blockRowPointers();
    const labelList& columnIndices = mbMatrix.blockColumnIndices();
    const scalarField& coeffs = mbMatrix.blockCoeffs();

    // Pout << "solve" << endl;
    // sleep(1);
    
    // Block size size
    MKL_INT blockSize = 6;

    // Number of cells
    label nCells = U.size();

    // 
    iparm40_ = mbMatrix.beam().csrAddr().globalNCellsOffset();
    iparm41_ = (mbMatrix.beam().csrAddr().globalNCellsOffset() + nCells) - 1;

    // Number of all cells
    globalNCells_ = sum(mbMatrix.beam().csrAddr().procNCells());

    // Matrix size (block matrix)
    MKL_INT n = globalNCells_;

    // RHS and solution vectors contruction
    // scalarField b(blockSize*n, 0);
    // scalarField x(blockSize*n, 0);
    scalarField& b = mbMatrix.rhs();
    scalarField& x = mbMatrix.solution();

    // // RHS and solution vector conversion
    // label k = 0;
    // for (label i=0; i<nCells; i++)
    // {
    //     for (label j=0; j<blockSize; j++)
    //     {
    //         b[k] = blockB[i](j);
    //         x[k] = U[i](j);
    //         k++;
    //     }
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
    int comm=0, rank=0, size=0;

    // Set number of threads
    mkl_set_num_threads(1);

    // Initialize MPI variables
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    comm = MPI_Comm_c2f(MPI_COMM_WORLD);

    // Setup Pardiso control parameters
    for ( label i=0; i<64; i++)
    {
        iparm[i] = 0;
    }
    iparm[ 0] =  1; /* Solver default parameters overriden with provided by iparm */
    iparm[ 1] =  2; /* Use METIS for fill-in reordering */
    iparm[ 5] =  0; /* Write solution into x */
    iparm[ 7] =  2; /* Max number of iterative refinement steps 2*/
    iparm[ 9] = 13; /* Perturb the pivot elements with 1E-13 */
    iparm[10] =  0; /* Use nonsymmetric permutation and scaling MPS */
    iparm[12] =  0; /* Switch on Maximum Weighted Matching algorithm (default for non-symmetric) */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[26] =  0; /* Check input data for correctness */
    iparm[34] =  1; /* Zero-based indexing of columns and rows */

    iparm[36] = blockSize; /* Block size */

    iparm[39] =  2; /* Input: matrix/rhs/solution are distributed between MPI processes  */
    /* If iparm[39]=2, the matrix is provided in distributed assembled matrix input          
       format. In this case, each MPI process stores only a part (or domain) of the matrix A 
       data. The bounds of the domain should be set via iparm(41) and iparm(42). Solution    
       vector is distributed between process in same manner with rhs. */
    iparm[40] = iparm40_; /* The number of row in global matrix, rhs element and solution vector
                      that begins the input domain belonging to this MPI process */
    iparm[41] = iparm41_; /* The number of row in global matrix, rhs element and solution vector
                      that ends the input domain belonging to this MPI process   */

    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum   = 1; /* Which factorization to use. */
    msglvl = 0; /* Print statistical information in file */
    error  = 0; /* Initialize error flag */

    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver
    for (label i=0; i<64; i++)
    {
        pt[i] = 0;
    }
    
    // Reordering and Symbolic Factorization. This step also allocates
    // all memory that is necessary for the factorization
    phase = 11;
    cluster_sparse_solver
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
        &comm,
        &error
    );

    // Pout << "solve 2" << endl;
    // sleep(1);
    
    // Pout << iparm[40] << ", " << iparm[41] << endl;
    // sleep(5);

    if (error!=0)
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector6>"
            "Foam::BlockMklClusterPardisoSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        ) << "ERROR during symbolic factorization: " << error
          << abort(FatalError);
    }

    // Numerical factorization
    phase = 22;
    cluster_sparse_solver
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
        &comm,
        &error
    );
    
    if (error != 0)
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector6>"
            "Foam::BlockMklClusterPardisoSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        ) << "ERROR during numerical factorization: " << error
          << abort(FatalError);
    }

    // Solve
    phase = 33;
    cluster_sparse_solver
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
        &comm,
        &error
    );
    
    if (error!=0)
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector6>"
            "Foam::BlockMklClusterPardisoSolver::solve"
            "("
            "    Field<Foam::vector6>& U,"
            "    const Field<Foam::vector6>& blockB"
            ")"
        ) << "ERROR during solution: " << error
          << abort(FatalError);
    }

    // Convert solution vector
    label k = 0;
    for (label i=0; i<nCells; i++)
    {
        for (label j=0; j<blockSize; j++)
        {
            U[i](j) = x[k++];
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

    // Termination and release of memory
    phase = -1; // Release internal memory
    cluster_sparse_solver
    (
        pt,
        &maxfct,
        &mnum,
        &mtype,
        &phase,
        &n,
        &ddum,
        &idum,
        &idum,
        &idum,
        &nrhs,
        iparm,
        &msglvl,
        &ddum,
        &ddum,
        &comm,
        &error
    );

    return solverPerf;
}


Foam::BlockSolverPerformance<Foam::vector6>
Foam::BlockMklClusterPardisoSolver::sequentialSolve
(
    Field<Foam::vector6>& U,
    const Field<Foam::vector6>& blockB
)
{
    // Prepare solver performance
    BlockSolverPerformance<vector6> solverPerf
    (
        typeName,
        this->fieldName()
    );

    vector6 norm = vector6::one; //this->normFactor(U, blockB);

    // labelList rowPointers;
    // labelList columnIndices;
    // scalarField coeffs;

    // convertFoamMatrixToMklClusterPardisoMatrix
    // (
    //     matrix_,
    //     rowPointers,
    //     columnIndices,
    //     coeffs
    // );

    multibeamFvBlockMatrix& mbMatrix =
        const_cast<multibeamFvBlockMatrix&>
        (
            static_cast<const multibeamFvBlockMatrix&>(matrix_)
        );

    const labelList& rowPointers = mbMatrix.blockRowPointers();
    const labelList& columnIndices = mbMatrix.blockColumnIndices();
    const scalarField& coeffs = mbMatrix.blockCoeffs();
    
    // Block size size
    MKL_INT blockSize = 6;

    // Number of cells
    MKL_INT nCells = U.size();

    // Matrix size (block matrix)
    MKL_INT n = nCells;

    // RHS and solution vectors contruction
    scalarField& b = mbMatrix.rhs();
    scalarField& x = mbMatrix.solution();
    
    // Calculate initial residual
    {
        scalarField p(x.size(), 0);
        MatVecMul
        (
            rowPointers,
            columnIndices,
            coeffs,
            x,
            p
        );

        Field<vector6> blockP(U.size());
        // matrix_.Amul(blockP, U);
        // Convert poroduct
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<blockSize; j++)
            {
                blockP[i](j) = p[k];
                k++;
            }
        }

        Field<vector6> blockR(blockB - blockP);

        solverPerf.initialResidual() =
            cmptDivide(gSum(cmptMag(blockR)), norm);
    }

    // scalarField b(blockSize*n, 0);
    // scalarField x(blockSize*n, 0);

    // // RHS and solution vector conversion
    // label k = 0;
    // for (label i=0; i<nCells; i++)
    // {
    //     for (label j=0; j<blockSize; j++)
    //     {
    //         b[k] = blockB[i](j);
    //         x[k] = U[i](j);
    //         k++;
    //     }
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
    iparm[1] = 2;          /* Fill-in reordering from METIS */
    iparm[3] = 0;          /* No iterative-direct algorithm */
    iparm[4] = 0;          /* No user fill-in reducing permutation */
    iparm[5] = 0;          /* Write solution into x */
    iparm[6] = 0;          /* Not in use */
    iparm[7] = 2;          /* Max numbers of iterative refinement steps 2 zt*/
    iparm[8] = 0;          /* Not in use */
    iparm[9] = 13;         /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 0;         /* Use nonsymmetric permutation and scaling MPS 1 zt*/
    iparm[11] = 0;         /* Conjugate transposed/transpose solve */
    iparm[12] = 0;         /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) 1 zt*/
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
    iparm[34] = 1;         /* Zero-based indexing of columns and rows */
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
    label k = 0;
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
        scalarField p(x.size(), 0);
        MatVecMul
        (
            rowPointers,
            columnIndices,
            coeffs,
            x,
            p
        );

        Field<vector6> blockP(U.size());
        // matrix_.Amul(blockP, U);
        // Convert poroduct
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<blockSize; j++)
            {
                blockP[i](j) = p[k];
                k++;
            }
        }

        Field<vector6> blockR(blockB - blockP);

        solverPerf.finalResidual() =
            cmptDivide(gSum(cmptMag(blockR)), norm);
        solverPerf.nIterations()++;
    }

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


void Foam::BlockMklClusterPardisoSolver::MatVecMul
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
