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
    BlockEigenSolverOF

Description
    Wrapper class to solve a 6x6 block coupled linear system using the Eigen
    header library.

    This class compiles with OpenFOAM-v2012.

Author
    Philip Cardiff, UCD
    Amirhossein Taran, UCD
    Colm McAlister UCD,
    Based on class from Zeljko Tukovic, FMENA, Zagreb

    
SourceFiles
    BlockEigenSolverOF.C

\*---------------------------------------------------------------------------*/

#ifndef BlockEigenSolverOF_H
#define BlockEigenSolverOF_H

#include "typeInfo.H"
#include "Switch.H"
#include "scalarMatrices.H"
#include "vector.H"
#include <Eigen/Sparse>
#include "RigidBodyCouplingData.H"


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

        //- Const-reference to owner list
        const labelList& own_;

        //- Const-reference to neighbour list
        const labelList& nei_;

        //- Optional monolithic rigid body coupling (rigidBodyEnd)
        RigidBodyMonolithicCoupling monolithicCoupling_;

        //- Newton increment of the body displacement from the last solve
        vector bodyTranslationIncrement_;

        //- Newton increment of the body step rotation vector from the last
        //  solve
        vector bodyRotationIncrement_;


    // Private Member Functions

        //- Number of 6x6 block rows: the beam cells, plus one for a
        //  monolithic rigid body
        label nBlocks() const
        {
            return d_.size() + (monolithicCoupling_.active ? 1 : 0);
        }

        //- Disallow default bitwise copy construct
        BlockEigenSolverOF(const BlockEigenSolverOF&);

        //- Disallow default bitwise assignment
        void operator=(const BlockEigenSolverOF&);


protected:

        //- Convert foam matrix to the Eigen matrix format
        void convertFoamMatrixToEigenMatrix
        (
            const Field<scalarSquareMatrix>& d,
            const Field<scalarSquareMatrix>& l,
            const Field<scalarSquareMatrix>& u,
            const labelList& own,
            const labelList& nei,
            Eigen::SparseMatrix<scalar>& A
        );


public:

    //- Runtime type information
    TypeName("Eigen");


    // Constructors

        //- Construct from matrix
        BlockEigenSolverOF
        (
            const Field<scalarSquareMatrix>& d,
            const Field<scalarSquareMatrix>& l,
            const Field<scalarSquareMatrix>& u,
            const labelList& own,
            const labelList& nei
        );


    // Destructor

        virtual ~BlockEigenSolverOF()
        {}


    // Member Functions

        //- Add a rigid body to the system, solved monolithically with the
        //  beam: six extra rows and columns for its translation and rotation
        //  Newton increments
        void setMonolithicCoupling(const RigidBodyMonolithicCoupling& coupling)
        {
            monolithicCoupling_ = coupling;
        }

        //- Solve the matrix with this source and solution vector and return the
        //  initial residual magnitude
        scalar solve
        (
            Field<scalarRectangularMatrix>& x,
            const Field<scalarRectangularMatrix>& b
        );

        //- Body displacement increment from the last solve (monolithic)
        const vector& bodyTranslationIncrement() const
        {
            return bodyTranslationIncrement_;
        }

        //- Body step rotation vector increment from the last solve
        //  (monolithic)
        const vector& bodyRotationIncrement() const
        {
            return bodyRotationIncrement_;
        }
};


} // End namespace Foam

#endif

// ************************************************************************* //
