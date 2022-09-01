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

#include "BlockPetscSolver.H"
#include "OFstream.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "multibeamFvBlockMatrix.H"
#include "beamModel.H"
#include "PstreamGlobals.H"


#include "petscSolver.H"
#include "petscControls.H"
#include "petscLinearSolverContexts.H"
#include "petscUtils.H"
#include "petscWrappedVector.H"
#include <cstring>

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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockPetscSolver::BlockPetscSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector6>& matrix,
    const dictionary& dict
)
:
    BlockIterativeSolver<vector6>(fieldName, matrix, dict),
    // BlockLduSolver<vector6>(fieldName, matrix, dict),
    petscDict_(dict.subDict("petsc")),
    eqName_(fieldName),
    prefix_("eqn_" + eqName_ + "_"),
    useBlockMatrix_(petscDict_.lookupOrDefault("useBlockMatrix", false)),
    iparm40_(-1),
    iparm41_(-1),
    globalNCells_(0),
    blockSize_(6)
{
    Info << "useBlockMatrix: "
         << useBlockMatrix_ << endl;
}


// ************************************************************************* //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector6>
Foam::BlockPetscSolver::solve
(
    Field<Foam::vector6>& U,
    const Field<Foam::vector6>& blockB
)
{
    // Prepare solver performance
    BlockSolverPerformance<vector6> blockSolverPerf
    (
        typeName,
        this->fieldName()
    );
    
    vector6 norm = vector6::one; //this->normFactor(U, blockB);

    multibeamFvBlockMatrix& mbMatrix =
        const_cast<multibeamFvBlockMatrix&>
        (
            static_cast<const multibeamFvBlockMatrix&>(matrix_)
        );

    if (useBlockMatrix_)
    {
        mbMatrix.rowMajor() = false;
    }
    
    const labelList& blockRowPointers = mbMatrix.blockRowPointers();
    const labelList& blockColumnIndices = mbMatrix.blockColumnIndices();
    const scalarField& blockCoeffs = mbMatrix.blockCoeffs();
    
    const labelList& rowPointers = mbMatrix.rowPointers();
    const labelList& columnIndices = mbMatrix.columnIndices();
    const scalarField& coeffs = mbMatrix.coeffs();

    // RHS and solution vectors contruction
    scalarField& source = mbMatrix.rhs();
    scalarField& psi = mbMatrix.solution();

    // Number of cells
    label nCells = U.size();


    // Petsc
    
    // Ensure PETSc is initialized
    const fvMesh& fvm = dynamicCast<const fvMesh>(matrix_.mesh().thisDb());
    const petscControls& controls = petscControls::New(fvm);

    if (!controls.valid())
    {
        FatalErrorInFunction
            << "PETSc not initialized" << nl << endl;
    }

    const bool verbose
    (
        petscDict_.lookupOrDefault("verbose", false)
    );

    // set all petsc flags in the prefix db
    dictionary petscDictOptions =
        petscDict_.subOrEmptyDict("options");
    PetscUtils::setFlags
    (
        prefix_,
        petscDictOptions,
        verbose
    );

    const petscLinearSolverContexts& contexts =
        petscLinearSolverContexts::New(fvm);

    petscLinearSolverContext& ctx = contexts.getContext(eqName_);
    const bool firsttimein = !ctx.initialized();

    dictionary petscDictCaching = petscDict_.subOrEmptyDict("caching");
    ctx.caching.init(petscDictCaching);
    ctx.caching.eventBegin();

    Mat& Amat = ctx.Amat;
    KSP& ksp  = ctx.ksp;

    // List<label>& lowNonZero = ctx.lowNonZero;
    // label& maxLowNonZeroPerRow = ctx.maxLowNonZeroPerRow;

    if (firsttimein)
    {
        DebugInfo<< "Initializing PETSc Linear Solver " << eqName_ << nl;

        ctx.initialized() = true;

        PetscLogStageRegister
        (
            ("foam_" + eqName_ + "_mat").c_str(),
            &ctx.matstage
        );
        PetscLogStageRegister
        (
            ("foam_" + eqName_ + "_pc").c_str(),
            &ctx.pcstage
        );
        PetscLogStageRegister
        (
            ("foam_" + eqName_ + "_ksp").c_str(),
            &ctx.kspstage
        );
        PetscLogStagePush(ctx.matstage);
        // buildMat(Amat, lowNonZero, maxLowNonZeroPerRow);

        // PetscBool flg;
        // MatGetOption(Amat, MAT_ROW_ORIENTED, &flg);
        // Info << "MAT_ROW_ORIENTED: " << flg << endl;

        if (useBlockMatrix_)
        {
            buildMatCSR
            (
                Amat,
                blockRowPointers,
                blockColumnIndices
            );
        }
        else
        {
            buildMatCSR
            (
                Amat,
                rowPointers,
                columnIndices
            );
        }

        PetscLogStagePop();
        buildKsp(Amat, ksp);

        // // No need for initial guess in the case of MUMPS direct solver
        // KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
    }

    const bool matup = ctx.caching.needsMatrixUpdate();
    const bool pcup = ctx.caching.needsPrecondUpdate();
    if (matup)
    {
        PetscLogStagePush(ctx.matstage);
        
        // updateMat(Amat, lowNonZero, maxLowNonZeroPerRow);

        if (useBlockMatrix_)
        {
            updateMatCSR
            (
                Amat,
                blockRowPointers,
                blockColumnIndices,
                blockCoeffs
            );
        }
        else
        {
            updateMatCSR
            (
                Amat,
                rowPointers,
                columnIndices,
                coeffs
            );
        }

        PetscLogStagePop();
    }
    updateKsp(ksp, Amat, !pcup);

    // Use built-in residual norm computation
    const bool usePetscResidualNorm
    (
        petscDict_.lookupOrDefault("use_petsc_residual_norm", false)
    );

    // Monitor built-in residual
    const bool monitorFoamResidualNorm
    (
        petscDict_.lookupOrDefault("monitor_foam_residual_norm", false)
    );

    // This optimization is disabled since some KSP implementations
    // are buggy in PETSc wrt KSP_NORM_NONE
    // users can still provide -ksp_norm_type none at command line
#if 0
    // Disable KSP default computation of residual norm if we are
    // monitoring convergence a-la OpenFOAM
    if (!usePetscResidualNorm)
    {
        KSPSetNormType(ksp, KSP_NORM_NONE);
    }
    else
    {
        KSPSetNormType(ksp, KSP_NORM_DEFAULT);
    }
#endif

    // ksp set options from db (may change norm type here if needed)
    if (firsttimein) KSPSetFromOptions(ksp);

    // Solver name from petsc (may have been changed from option database)
    KSPType ksptype;
    KSPGetType(ksp, &ksptype);
    word solverName(ksptype);

    // if (solverName == word("preonly"))
    // {
    //     PC pc;
    //     KSPGetPC(ksp, &pc);
    //     const char* type;
    //     PCFactorGetMatSolverType(pc, &type);
    //     solverName += word(type);
    // }

    if (firsttimein)
    {
        if (solverName == word("preonly"))
        {
            // No need for initial guess in the case of MUMPS direct solver
            KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
        }
    }
    
    // Setup class containing solver performance data
    lduSolverPerformance solverPerf
    (
        "PETSc-" + solverName,
        fieldName()
    );

    blockSolverPerf.solverName() = solverPerf.solverName();
    
    // Retain copy of solverPerformance
    ctx.performance = solverPerf;

    // Create solution and rhs vectors for PETSc
    // We could create these once, and call PlaceArray/ResetArray instead
    PetscWrappedVector petsc_psi(psi, Amat);
    PetscWrappedVector petsc_source(source, Amat);

    // Initialize data to compute L1 residual or to monitor it while solving
    ctx.performance.initialResidual() = 0;
    KSPNormType normType;
    KSPGetNormType(ksp, &normType);
    if
    (
       !usePetscResidualNorm
     || monitorFoamResidualNorm
     || normType == KSP_NORM_NONE
    )
    {
        ctx.createAuxVecs();
        ctx.initArowsSumVec();
        ctx.computeNormFactor(petsc_psi,petsc_source);
        if (normType == KSP_NORM_NONE && !usePetscResidualNorm)
        {
           ctx.performance.initialResidual() =
               ctx.getResidualNormL1(petsc_source);
        }
    }

    // Convergence testing
    if (firsttimein)
    {
        if (usePetscResidualNorm)
        {
            // Add monitor to record initial residual computed by PETSc
            KSPMonitorSet
            (
                ksp,
                PetscUtils::foamKSPMonitorRecordInit,
                &ctx,
                NULL
            );
        }
        else
        {
            ctx.useFoamTest = true;
            KSPSetConvergenceTest
            (
                ksp,
                PetscUtils::foamKSPConverge,
                &ctx,
                NULL
            );
        }

        // Monitoring (print to stdout) residual reduction in OpenFOAM norm
        if (monitorFoamResidualNorm)
        {
            KSPMonitorSet
            (
                ksp,
                PetscUtils::foamKSPMonitorFoam,
                &ctx,
                NULL
            );
        }
    }

    // Geometric multigrid support
    // (pc_type ml or pc_type gamg pc_gamg_type geo for 2d only as for 3.14)
    if (firsttimein)
    {
        PC pc;
        void (*f)(void) = NULL;

        KSPGetPC(ksp, &pc);
        PetscObjectQueryFunction((PetscObject)pc, "PCSetCoordinates_C", &f);

        // Cell-centres are contiguous in memory (array of structures)
        // Not sure about const-ness here
        // Could also use the PrecisionAdaptor
        if (f)
        {
            const vectorField& cc = fvm.C().internalField();
            // const vectorField& cc = fvm.C().primitiveField();
            const PetscInt n = cc.size();
            const PetscInt sdim = vector::nComponents;

            List<PetscReal> ccPoints(cc.size()*vector::nComponents);

            auto iter = ccPoints.data();
            for (const vector& v : cc)
            {
                *(iter++) = v.x();
                *(iter++) = v.y();
                *(iter++) = v.z();
            }

            PCSetCoordinates(pc, sdim, n, ccPoints.data());
        }
    }

    // Setup KSP (this is not explicitly needed, but we call it to separate
    // PCSetUp from KSPSolve timings when requesting -log_view from PETSc)
    PetscLogStagePush(ctx.pcstage);
    KSPSetUp(ksp);
    KSPSetUpOnBlocks(ksp);
    PetscLogStagePop();

    // Calculate initial residual
    {
        Vec res;
        MatCreateVecs(Amat, NULL, &res);
        MatMult(Amat, petsc_psi, res);
        VecAXPY(res, -1., petsc_source);

        PetscScalar *array;
        VecGetArray(res, &array);

        Field<vector6> blockRes(U.size());
        
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<blockSize_; j++)
            {
                blockRes[i](j) = array[k++];
            }
        }
        
        VecRestoreArray(res, &array);

        blockSolverPerf.initialResidual() =
            cmptDivide(gSum(cmptMag(blockRes)), norm);
    }

    // If using the PETSc infrastructure to monitor convergence,
    // we need to report back the initial residual
    if (usePetscResidualNorm)
    {
        PetscReal rnorm;
        KSPNormType normType;
        KSPGetNormType(ksp, &normType);
        if (normType == KSP_NORM_NONE)
        // We run with PETSc norm KSP_NORM_NONE -> report L1 norm
        {
            rnorm = ctx.getResidualNormL1(petsc_psi, petsc_source);
        }
        else
        // report final PETSc norm since users explicitly ask
        // to run with PETSc norm
        {
           KSPGetResidualNorm(ksp, &rnorm);
        }
        ctx.performance.finalResidual() = rnorm;

        // Info << "irnorm: " << rnorm << endl;
    }
    
    // Solve A x = b
    PetscLogStagePush(ctx.kspstage);
    KSPSolve(ksp, petsc_source, petsc_psi);
    PetscLogStagePop();

    ctx.caching.eventEnd();

    // Set nIterations and final residual
    PetscInt nIters;
    KSPGetIterationNumber(ksp, &nIters);
    ctx.performance.nIterations() = nIters;
    // Info << "nIters: " << nIters << endl;

    // If using the PETSc infrastructure to monitor convergence,
    // we need to report back the final residual
    if (usePetscResidualNorm)
    {
        PetscReal rnorm;
        KSPNormType normType;
        KSPGetNormType(ksp, &normType);
        if (normType == KSP_NORM_NONE)
        // We run with PETSc norm KSP_NORM_NONE -> report L1 norm
        {
            rnorm = ctx.getResidualNormL1(petsc_psi, petsc_source);
        }
        else
        // report final PETSc norm since users explicitly ask
        // to run with PETSc norm
        {
           KSPGetResidualNorm(ksp, &rnorm);
        }
        ctx.performance.finalResidual() = rnorm;

        // Info << "rnorm: " << rnorm << endl;
    }

    // Convert solution vector
    label k = 0;
    for (label i=0; i<nCells; i++)
    {
        for (label j=0; j<blockSize_; j++)
        {
            U[i](j) = psi[k++];
        }
    }

    // Calculate final residual
    {
        Vec res;
        MatCreateVecs(Amat, NULL, &res);
        MatMult(Amat, petsc_psi, res);
        VecAXPY(res, -1., petsc_source);

        PetscScalar *array;
        VecGetArray(res, &array);

        Field<vector6> blockRes(U.size());
        
        label k = 0;
        for (label i=0; i<nCells; i++)
        {
            for (label j=0; j<blockSize_; j++)
            {
                blockRes[i](j) = array[k++];
            }
        }

        VecRestoreArray(res, &array); 

        blockSolverPerf.finalResidual() =
            cmptDivide(gSum(cmptMag(blockRes)), norm);

        blockSolverPerf.nIterations() = nIters;
    }

    return blockSolverPerf;
}


void Foam::BlockPetscSolver::updateKsp
(
    KSP& ksp,
    Mat& Amat,
    const bool reuse
) const
{
    if (reuse)
    {
        DebugInfo<< "Cache-Hit: reuse preconditioner " << eqName_ << nl;

        KSPSetReusePreconditioner(ksp, PETSC_TRUE);
    }
    else
    {
        DebugInfo<< "Cache-Hit: rebuild preconditioner " << eqName_ << nl;

        KSPSetReusePreconditioner(ksp, PETSC_FALSE);
    }

    KSPSetOperators(ksp, Amat, Amat);

    // update tolerance and relTol for the (*)Final solver
    // Info << relTolerance() << ", " << tolerance() << ", "
    //      << maxIter() << endl;    
    KSPSetTolerances
    (
        ksp,
        relTolerance(),
        tolerance(),
        PETSC_DEFAULT,
        maxIter()
    );
}


void Foam::BlockPetscSolver::buildKsp
(
    Mat& Amat,
    KSP& ksp
) const
{
    // Create parallel solver context
    KSPCreate(PetscObjectComm((PetscObject)Amat), &ksp);

    // Set the prefix for the options db (e.g. -eqn_p_)
    KSPSetOptionsPrefix(ksp, prefix_.c_str());

    // ksp set operator and preconditioner
    KSPSetOperators(ksp, Amat, Amat);

    // OpenFOAM relative tolerance -> tolerance_
    KSPSetTolerances
    (
        ksp,
        relTolerance(),
        tolerance(),
        PETSC_DEFAULT,
        maxIter()
    );

    // Use solution from the previous timestep as initial guess
    KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
}


void Foam::BlockPetscSolver::buildMatCSR
(
    Mat& Amat,
    const labelList& rowPointers,
    const labelList& columnIndices,
    const bool symmetric
) const
{
    // communicator for the matrix
    // PETSc internally duplicates the communicator to
    // avoid tag/attributes clashing

    const auto comm =
    (
        PstreamGlobals::MPICommunicators_.empty()
      ? PETSC_COMM_SELF
      : PstreamGlobals::MPICommunicators_[matrix_.mesh().comm()]
    );

    // Local degrees-of-freedom i.e. number of local rows
    label nRows = (rowPointers.size()-1);
    if (useBlockMatrix_)
    {
        nRows *= blockSize_;
    }

    // Create a petsc matrix
    MatCreate(comm, &Amat);
    MatSetSizes(Amat, nRows, nRows, PETSC_DETERMINE, PETSC_DETERMINE);

    // Set the prefix for the options db (e.g. -eqn_p_)
    MatSetOptionsPrefix(Amat, prefix_.c_str());
    MatSetFromOptions(Amat);

    // Preallocate the matrix using CSR addressing
    
    if (useBlockMatrix_)
    {
        MatMPIBAIJSetPreallocationCSR
        (
            Amat,
            blockSize_,
            rowPointers.cdata(),
            columnIndices.cdata(),
            PETSC_NULL
        );
    }
    else
    {
        MatMPIAIJSetPreallocationCSR
        (
            Amat,
            rowPointers.cdata(),
            columnIndices.cdata(),
            PETSC_NULL
        );
    }

    // set the matrix options
    if (symmetric)
    {
        MatSetOption(Amat, MAT_SYMMETRIC, PETSC_TRUE);
        MatSetOption(Amat, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);
    }
    else
    {
        MatSetOption(Amat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);
    }

    MatSetOption(Amat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);

    // Skip MPI_Allreduce calls when finalizing assembly
    MatSetOption(Amat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);

    // Allow general setup of other matrices instead of erroring
    MatSetUp(Amat);
}


void Foam::BlockPetscSolver::updateMatCSR
(
    Mat& Amat,
    const labelList& rowPointers,
    const labelList& columnIndices,
    const scalarField& coeffs
) const
{
    const auto comm =
    (
        PstreamGlobals::MPICommunicators_.empty()
      ? PETSC_COMM_SELF
      : PstreamGlobals::MPICommunicators_[matrix_.mesh().comm()]
    );

    // Local degrees-of-freedom i.e. number of local rows
    label nRows = (rowPointers.size()-1);
    if (useBlockMatrix_)
    {
        nRows *= blockSize_;
    }

    if (useBlockMatrix_)
    {
        MatCreateMPIBAIJWithArrays
        (
            comm,
            blockSize_,
            nRows,
            nRows,
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            rowPointers.cdata(),
            columnIndices.cdata(),
            coeffs.cdata(),
            &Amat
        );
    }
    else
    {
        MatCreateMPIAIJWithArrays
        (
            comm,
            nRows,
            nRows,
            PETSC_DETERMINE,
            PETSC_DETERMINE,
            rowPointers.cdata(),
            columnIndices.cdata(),
            coeffs.cdata(),
            &Amat
        );
    }

    MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY);
}
