/*---------------------------------------------------------------------------*\
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
    Foam::fixedDisplacementFvPatchVectorField

Description
    Fixed displacement boundary condition for beams

SourceFiles
    fixedDisplacementFvPatchVectorField.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementFvPatchVectorField_H
#define fixedDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedDisplacementFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Constant total displacement field
        vectorField totalDisp_;

        //- Displacement time series
        interpolationTable<vector> dispSeries_;

    // Private Functions

public:

    //- Runtime type information
    TypeName("fixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedDisplacementFvPatchVectorField
        //  onto a new patch
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedDisplacementFvPatchVectorField
        (
            const fixedDisplacementFvPatchVectorField&,
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
                new fixedDisplacementFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access functions

        //- Update the coefficients associated with the patch field
            virtual vectorField& totalDisp()
            {
                return totalDisp_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
