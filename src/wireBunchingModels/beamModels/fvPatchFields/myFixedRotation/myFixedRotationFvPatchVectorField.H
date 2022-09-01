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
    myFixedRotationFvPatchVectorField

Description
    Specifying a beam rotation value at the patch

SourceFiles
    myFixedRotationFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef myFixedRotationFvPatchVectorField_H
#define myFixedRotationFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class myFixedRotationFvPatch Declaration
\*---------------------------------------------------------------------------*/

class myFixedRotationFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data
        
        //- Constant total rotation field
        vectorField totalTheta_;

        //- Theta time series
        interpolationTable<vector> thetaSeries_;
        
        //- Current time index
        label timeIndex_;

public:

    //- Runtime type information
    TypeName("myFixedRotation");


    // Constructors

        //- Construct from patch and internal field
        myFixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myFixedRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  myFixedRotationFvPatchVectorField onto a new patch
        myFixedRotationFvPatchVectorField
        (
            const myFixedRotationFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myFixedRotationFvPatchVectorField
        (
            const myFixedRotationFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new myFixedRotationFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myFixedRotationFvPatchVectorField
        (
            const myFixedRotationFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new myFixedRotationFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        //- Update the coefficients associated with the patch field
            virtual vectorField& totalTheta()
            {
                return totalTheta_;
            }
            
        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // //- Evaluate the patch field
        // virtual void evaluate
        // (
        //     const Pstream::commsTypes commsType = Pstream::blocking
        // );
  
        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
