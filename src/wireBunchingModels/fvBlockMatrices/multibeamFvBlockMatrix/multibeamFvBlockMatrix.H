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
    multibeamFvBlockMatrix

Description

Author
    Zeljko Tukovic, FMENA, Zagreb 

SourceFiles
    multibeamFvBlockMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef multibeamFvBlockMatrix_H
#define multibeamFvBlockMatrix_H

#include "fvc.H"
#include "fvBlockMatrix.H"
#include "typeInfo.H"
#include "Switch.H"
#include "vector6.H"
#include "tensor6.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class beamModel;

/*---------------------------------------------------------------------------* \
                      Class multibeamFvBlockMatrix Declaration
\*---------------------------------------------------------------------------*/

typedef Pair<tensor6> tensor6Pair;

typedef List<tensor6Pair> tensor6PairList;
typedef List<tensor6PairList> tensor6PairListList;

typedef fvBlockMatrix<vector6> fvBlockVector6Matrix;

class multibeamFvBlockMatrix
:
    public fvBlockVector6Matrix
{
    // Private data

        //- All beams in one model
        beamModel& beam_;

        //- Line contact coefficients
        mutable PtrList<tensor6PairListList> lineContactCoeffs_;

        //- Point contact coefficients
        mutable tensor6PairListList pointContactCoeffs_;

        //- Active conical pulley contacts
        Map<FixedList<label, 3> > activeConicalPulleyContacts_;

        //- CSR addressing row pointers
        mutable labelList* blockRowPointersPtr_;

        //- CSR addressing column indices
        mutable labelList* blockColumnIndicesPtr_;

        //- CSR matrix format coefficients
        mutable scalarField* blockCoeffsPtr_;

        //- CSR addressing row pointers
        mutable labelList* rowPointersPtr_;

        //- CSR addressing column indices
        mutable labelList* columnIndicesPtr_;

        //- CSR matrix format coefficients
        mutable scalarField* coeffsPtr_;

        //- CSR matrix format right hand side vector
        mutable scalarField* rhsPtr_;

        //- CSR matrix format solution vector
        mutable scalarField* solutionPtr_;

        //- blockCoeffs filling order
        mutable bool rowMajor_;

    // Private Member Functions

        //- Convert LDU (OF) to CSR matrix
        void makeBlockCsrMatrix() const;

        //- Convert LDU (OF) to CSR matrix
        void makeCsrMatrix() const;
    
        //- Find point contact
        label upperPointContact(const label bI, const label segI) const;
    
        //- Find point contact
        label lowerPointContact(const label bI, const label segI) const;
    
        //- Disallow default bitwise copy construct
        multibeamFvBlockMatrix(const multibeamFvBlockMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const multibeamFvBlockMatrix&);

protected:

public:

    //- Runtime type information
    TypeName("multibeamFvBlockMatrix");

    // Constructors

        //- Construct given a field to solve for
        multibeamFvBlockMatrix
        (
            GeometricField<vector6, fvPatchField, volMesh>& phi,
            beamModel& beam
        );

        // //- Construct as copy
        // multibeamFvBlockMatrix(const fvBlockMatrix<vector6>&);


    //- Destructor
    virtual ~multibeamFvBlockMatrix();

    
    // Member Functions

        label nCells() const
        {
            return this->lduAddr().size();
        }

        label blockSize() const
        {
            return 6;
        }
    
        bool hybrid() const
        {
            return (activeConicalPulleyContacts().size() > 0);
        }
    
        const beamModel& beam() const
        {
            return beam_;
        }
    
        beamModel& beam()
        {
            return beam_;
        }
    
        //- Return line contact off-diagonal coefficients
        const PtrList<tensor6PairListList>& lineContactCoeffs() const
        {
            return lineContactCoeffs_;
        }

        //- Return line contact off-diagonal coefficients
        PtrList<tensor6PairListList>& lineContactCoeffs();

        //- Return point contact off-diagonal coefficients
        const tensor6PairListList& pointContactCoeffs() const
        {
            return pointContactCoeffs_;
        }

        //- Return point contact off-diagonal coefficients
        tensor6PairListList& pointContactCoeffs();

        //- Return active conical pulley contacts
        Map<FixedList<label, 3> >& activeConicalPulleyContacts()
        {
            return activeConicalPulleyContacts_;
        }

        //- Return active conical pulley contacts
        const Map<FixedList<label, 3> >&
        activeConicalPulleyContacts() const
        {
            return activeConicalPulleyContacts_;
        }

        //- Return CSR matrix row pointers
        const labelList& blockRowPointers() const;

        //- Return CSR matrix column indices
        const labelList& blockColumnIndices() const;

        //- Return CSR matrix coefficients
        const scalarField& blockCoeffs() const;

        //- Return CSR matrix row pointers
        const labelList& rowPointers() const;

        //- Return CSR matrix column indices
        const labelList& columnIndices() const;

        //- Return CSR matrix coefficients
        const scalarField& coeffs() const;
    
        //- Return CSR matrix right hand side vector
        const scalarField& rhs() const;
    
        //- Return CSR matrix right hand side vector
        scalarField& rhs();
    
        //- Return CSR matrix solution vector
        const scalarField& solution() const;
    
        //- Return CSR matrix solution vector
        scalarField& solution();

        //- Retrive lagrange multipliers from the solution vector
        void retriveLagrangeMultipliers();

        //- blockCoeffs filling order
        bool rowMajor() const
        {
            return rowMajor_;
        }
    
        //- blockCoeffs filling order
        bool& rowMajor()
        {
            Info << "blockRowPointersPtr_: " << blockRowPointersPtr_ << endl;
            return rowMajor_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
