/*---------------------------------------------------------------------------*\
    Copyright (C) Hrvoje Jasak, Wikki Ltd
-------------------------------------------------------------------------------
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    BlockTDMSolver

Description
    Preconditioned Bi-Conjugate Gradient stabilised solver.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    BlockTDMSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockTDMSolver_H
#define BlockTDMSolver_H

#include "BlockBiCGStabSolver.H"
#include "vector6.H"
#include "tensor6.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class BlockTDMSolver Declaration
\*---------------------------------------------------------------------------*/

typedef BlockBiCGStabSolver<vector6> blockBiCGStabVector6Solver;

class BlockTDMSolver
:
    public BlockIterativeSolver<vector6>
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockTDMSolver(const BlockTDMSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockTDMSolver&);


public:

    //- Runtime type information
    TypeName("BlockTDM");


    // Constructors

        //- Construct from matrix components and solver data stream
        BlockTDMSolver
        (
            const word& fieldName,
            const BlockLduMatrix<vector6>& matrix,
            const dictionary& dict
        );


        // Destructor

            virtual ~BlockTDMSolver()
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
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
