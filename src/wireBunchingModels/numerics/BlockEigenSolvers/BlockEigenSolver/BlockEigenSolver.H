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
    BlockEigenSolver

Description

Author
    Zeljko Tukovic, FMENA, Zagreb 

SourceFiles
    BlockEigenSolver.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenSolver_H
#define BlockEigenSolver_H

#include <functional>
#include <numeric>

#include "typeInfo.H"
#include "Switch.H"
#include "vector6.H"
#include "tensor6.H"
#include "BlockLduSolver.H"

#include <Eigen/Sparse>
// #include "BlockSparseMatrix.h"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockEigenSolver Declaration
\*---------------------------------------------------------------------------*/

typedef BlockLduSolver<vector6> blockVector6Solver;

class BlockEigenSolver
:
    public blockVector6Solver
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockEigenSolver(const BlockEigenSolver&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenSolver&);

protected:

        //- Convert foam matrix to the MKL-PARDISO matrix format
        void convertFoamMatrixToEigenMatrix
        (
            const BlockLduMatrix<vector6>& matrix,
            Eigen::SparseMatrix<scalar>& A,
            Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b,
            Eigen::Matrix<scalar, Eigen::Dynamic, 1>& x
        );

public:

    //- Runtime type information
    TypeName("Eigen");

    // Constructors

        //- Construct from matrix
        BlockEigenSolver
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
