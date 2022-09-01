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

Class
    BlockPetscSolver

Description

Author
    Zeljko Tukovic, FMENA, Zagreb 

SourceFiles
    BlockPetscSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockPetscSolver_H
#define BlockPetscSolver_H

#include "typeInfo.H"
#include "Switch.H"
#include "vector6.H"
#include "tensor6.H"
#include "BlockLduSolver.H"
#include "BlockIterativeSolver.H"

// PETSc
#include "petscmat.h"
#include "petscksp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
/*---------------------------------------------------------------------------*\
                      Class BlockPetscSolver Declaration
\*---------------------------------------------------------------------------*/

typedef BlockIterativeSolver<vector6> blockVector6Solver;
// typedef BlockLduSolver<vector6> blockVector6Solver;
  
class BlockPetscSolver
:
    public blockVector6Solver
{
    // Private data

        //- Coefficients for the petsc solvers
        dictionary petscDict_;

        //- Equation name for solver context
        word eqName_;

        //- Solution name used as prefix
        word prefix_;

        //- Use block matrix
        bool useBlockMatrix_;
        
        //- Global index of first equation on this processor
        label iparm40_;
    
        //- Global index of last equation on this processor
        label iparm41_;

        //- Global number of cells
        label globalNCells_;

        //- Block size
        label blockSize_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockPetscSolver(const BlockPetscSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockPetscSolver&);

        //- Build the ksp_ krylov solver context
        void buildKsp(Mat& Amat, KSP& ksp) const;

        //- Update the ksp_ krylov solver context
        void updateKsp(KSP& ksp, Mat& Amat, const bool reuse) const;

        //- Compute allocation (proc and off-proc)
        //- for the petsc matrix Amat
        void buildMatCSR
        (
            Mat& Amat,
            const labelList& rowPointers,
            const labelList& columnIndices,
            const bool symmetric = false
        ) const;

        //- Update Amat matrix by inserting new values
        void updateMatCSR
        (
            Mat& Amat,
            const labelList& rowPointers,
            const labelList& columnIndices,
            const scalarField& coeffs
        ) const;

protected:

        // //- Convert foam matrix to the MKL-PARDISO matrix format
        // void convertFoamMatrixToPetscMatrix
        // (
        //     const BlockLduMatrix<vector6>& matrix,
        //     labelList& rowPointers,
        //     labelList& columnIndices,
        //     scalarField& coeffs
        // );
    
public:

    //- Runtime type information
    TypeName("BlockPetsc");

    // Constructors

        //- Construct from matrix
        BlockPetscSolver
        (
            const word& fieldName,
            const BlockLduMatrix<vector6>& matrix,
            const dictionary& dict
        );

    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<vector6> solve
        (
            Field<vector6>& x,
            const Field<vector6>& b
        );

        void MatVecMul
        (
            const labelList& rowPointers,
            const labelList& columnIndices,
            const scalarField& coeffs,
            const scalarField& x,
            scalarField& p
        ) const;

        //- Re-initialise solver after matrix coefficient update
        virtual void initMatrix()
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
