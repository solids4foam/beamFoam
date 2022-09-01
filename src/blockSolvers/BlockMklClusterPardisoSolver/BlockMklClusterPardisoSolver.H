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
    BlockMklClusterPardisoSolver

Description

Author
    Zeljko Tukovic, FMENA, Zagreb 

SourceFiles
    BlockMklClusterPardisoSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockMklClusterPardisoSolver_H
#define BlockMklClusterPardisoSolver_H

#include "typeInfo.H"
#include "Switch.H"
#include "vector6.H"
#include "tensor6.H"
#include "BlockLduSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class BlockMklClusterPardisoSolver Declaration
\*---------------------------------------------------------------------------*/

typedef BlockLduSolver<vector6> blockVector6Solver;
  
class BlockMklClusterPardisoSolver
:
    public blockVector6Solver
{
    // Private data

       //- Global index of first equation on this processor
       label iparm40_;
    
       //- Global index of last equation on this processor
       label iparm41_;

       //- Global number of cells
       label globalNCells_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockMklClusterPardisoSolver(const BlockMklClusterPardisoSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockMklClusterPardisoSolver&);

protected:

        //- Convert foam matrix to the MKL-PARDISO matrix format
        void convertFoamMatrixToMklClusterPardisoMatrix
        (
            const BlockLduMatrix<vector6>& matrix,
            labelList& rowPointers,
            labelList& columnIndices,
            scalarField& coeffs
        );

public:

    //- Runtime type information
    TypeName("MklClusterPardiso");

    // Constructors

        //- Construct from matrix
        BlockMklClusterPardisoSolver
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

        //- Solve the matrix with this solver
        BlockSolverPerformance<vector6> sequentialSolve
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
