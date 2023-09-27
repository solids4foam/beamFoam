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
    BlockEigenSolverOF

Description
    Wrapper class to solve a 6x6 block coupled linear system using the Eigen
    header library.

    This class compiles with OpenFOAM-v2012.

Author
    Philip Cardiff, UCD
    Amirhossein Taran, UCD
    Based on class from Zeljko Tukovic, FMENA, Zagreb 

SourceFiles
    BlockEigenSolverOF.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenSolverOF_H
#define BlockEigenSolverOF_H

// #include <functional>
// #include <numeric>
#include "typeInfo.H"
#include "Switch.H"
#include <Eigen/Sparse>
#include "scalarMatrices.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BlockEigenSolverOF Declaration
\*---------------------------------------------------------------------------*/

class BlockEigenSolverOF
{
    // Private Data

        //- Const-reference to the sparse matrix diagonal
        const Field<scalarSquareMatrix>& d_;

        //- Const-reference to the sparse matrix upper
        const Field<scalarSquareMatrix>& l_;

        //- Const-reference to the sparse matrix lower
        const Field<scalarSquareMatrix>& u_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        BlockEigenSolverOF(const BlockEigenSolverOF&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenSolverOF&);

protected:

        //- Convert foam matrix to the MKL-PARDISO matrix format
        // void convertFoamMatrixToEigenMatrix
        // (
        //     const BlockLduMatrix<vector6>& matrix,
        //     Eigen::SparseMatrix<scalar>& A,
        //     Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b,
        //     Eigen::Matrix<scalar, Eigen::Dynamic, 1>& x
        // );

public:

    //- Runtime type information
    TypeName("Eigen");

    // Constructors

        //- Construct from matrix
        BlockEigenSolverOF
        (
            const Field<scalarSquareMatrix>& d,
            const Field<scalarSquareMatrix>& l,
            const Field<scalarSquareMatrix>& u
            // const dictionary& dict
        );

    // Destructor

        virtual ~BlockEigenSolverOF()
        {}

    // Member Functions

        //- Solve the matrix with this solver
        virtual bool solve
        (
            Field<scalarRectangularMatrix>& x,
            const Field<scalarRectangularMatrix>& b
        );

        // //- Re-initialise solver after matrix coefficient update
        // virtual void initMatrix()
        // {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
