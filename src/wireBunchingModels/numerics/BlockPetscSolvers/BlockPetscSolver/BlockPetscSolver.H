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
    An abstract base class wrapper for linear solvers from the Eigen
    (http://eigen.tuxfamily.org) library.

    Notes
    - The solver is currently only implemented for vector6
    - I would recommend removing '-Wold-style-cast' from
      $WM_PROJECT_DIR/wmake/rules/$WM_ARCH/c++ as the terminal will be polluted
      with 'old-stlye-casts' warnings from the compilation of the Eigen header
      files.

Author
    Zeljko Tukovic, FSB Zagreb
    Philip Cardiff, UCD.  All rights reserved

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

#include "petsc.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class BlockPetscSolver Declaration
\*---------------------------------------------------------------------------*/

typedef BlockLduSolver<vector6> blockVector6Solver;
  
class BlockPetscSolver
:
    public blockVector6Solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockPetscSolver(const BlockPetscSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockPetscSolver&);

protected:

        //- Check if the matrix is for a 2-D case
        bool checkTwoD() const;

        //- Calculate the number of degrees of freedom
        int calcDegreesOfFreedom
        (
            const BlockLduMatrix<vector6>& matrix,
            const bool twoD
        ) const;

        //- Convert foam matrix to the Eigen matrix format
        void convertFoamMatrixToPetscMatrix
        (
            const BlockLduMatrix<vector6>& matrix,
            Mat& A
            // Eigen::SparseMatrix<scalar>& A
        );

public:

    //- Runtime type information
    TypeName("PetscSolver");

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

        //- Re-initialise solver after matrix coefficient update
        virtual void initMatrix()
        {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
