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

Class
    implBeamsContactBlockEigenSparseLUSolver

Description
    The SparseLU direct sparse linear equation solver from package Eigen is
    called to solve the system.

    Notes
    - The solver does not run in parallel
    - The solver is currently only implemented for two vectors
    - I would recommend removing '-Wold-style-cast' from
      $WM_PROJECT_DIR/wmake/rules/$WM_ARCH/c++ as the terminal will be polluted
      with 'old-stlye-casts' warnings from the compilation of the Eigen header
      files.

Author
    Philip Cardiff, UCD.  All rights reserved

SourceFiles
    implBeamsContactBlockEigenSparseLUSolver.C

\*---------------------------------------------------------------------------*/

#ifndef implBeamsContactBlockEigenSparseLUSolver_H
#define implBeamsContactBlockEigenSparseLUSolver_H

#include "BlockEigenSolver.H"

#include "Field.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"
#include "BlockLduMatrix.H"
#include "BlockLduSolver.H"
#include "blockLduSolvers.H"
#include "BlockSolverPerformance.H"
#include <Eigen/Sparse>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
         Class implBeamsContactBlockEigenSparseLUSolver Declaration
\*---------------------------------------------------------------------------*/

class implBeamsContactBlockEigenSparseLUSolver
:
    public BlockEigenSolver
{
  
  
    // Private Member Functions

        //- Disallow default bitwise copy construct
        implBeamsContactBlockEigenSparseLUSolver
        (
            const implBeamsContactBlockEigenSparseLUSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const implBeamsContactBlockEigenSparseLUSolver&);

public:

    //- Runtime type information
    TypeName("implBeamsContactEigenSparseLU");

    // Constructors

        //- Construct from matrix
        implBeamsContactBlockEigenSparseLUSolver
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

        //- Add implicit contact to the matrix
        void addImplicitContact
        (
            Eigen::SparseMatrix<scalar>& A,
            Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b
        );
  
        // //- Return contact force as a function of gap 
        // scalar contactForce(const scalar g) const;
  
        // //- Return contact force derivative as a function of gap 
        // scalar contactForceDerivative(const scalar g) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
