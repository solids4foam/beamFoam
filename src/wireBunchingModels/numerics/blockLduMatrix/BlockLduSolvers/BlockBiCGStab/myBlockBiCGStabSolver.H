/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    myBlockBiCGStabSolver

Description
    Preconditioned Bi-Conjugate Gradient stabilised solver.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    myBlockBiCGStabSolver.C

\*---------------------------------------------------------------------------*/

#ifndef myBlockBiCGStabSolver_H
#define myBlockBiCGStabSolver_H

#include "BlockBiCGStabSolver.H"
#include "vector6.H"
#include "tensor6.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class myBlockBiCGStabSolver Declaration
\*---------------------------------------------------------------------------*/

typedef BlockBiCGStabSolver<vector6> blockBiCGStabVector6Solver;

class myBlockBiCGStabSolver
:
    public blockBiCGStabVector6Solver
{
    // Private data

        //- Preconditioner
        autoPtr<BlockLduPrecon<vector6> > preconPtr_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        myBlockBiCGStabSolver(const myBlockBiCGStabSolver&);

        //- Disallow default bitwise assignment
        void operator=(const myBlockBiCGStabSolver&);


public:

    //- Runtime type information
    TypeName("myBiCGStab");


    // Constructors

        //- Construct from matrix components and solver data stream
        myBlockBiCGStabSolver
        (
            const word& fieldName,
            const BlockLduMatrix<vector6>& matrix,
            const dictionary& dict
        );


        // Destructor

            virtual ~myBlockBiCGStabSolver()
            {}


    // Member Functions

        //- Solve the matrix with this solver
        virtual BlockSolverPerformance<vector6> solve
        (
            Field<vector6>& x,
            const Field<vector6>& b
        );

        //- Re-initialise solver after matrix coefficient update
        virtual void initMatrix()
        {
            preconPtr_->initMatrix();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
