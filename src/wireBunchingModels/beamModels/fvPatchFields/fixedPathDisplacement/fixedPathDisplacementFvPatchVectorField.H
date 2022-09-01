/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::fixedPathDisplacementFvPatchVectorField

Description
    Fixed displacement boundary condition for beams

SourceFiles
    fixedPathDisplacementFvPatchVectorField.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef fixedPathDisplacementFvPatchVectorField_H
#define fixedPathDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "HermiteSpline.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fixedPathDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class fixedPathDisplacementFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- File name
        fileName fileName_;

        //- Path defined by Hermite spline
        autoPtr<HermiteSpline> pathPtr_;

        //- Velocity
        scalar velocity_;

    // Private Functions

public:

    //- Runtime type information
    TypeName("fixedPathDisplacement");


    // Constructors

        //- Construct from patch and internal field
        fixedPathDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedPathDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedPathDisplacementFvPatchVectorField
        //  onto a new patch
        fixedPathDisplacementFvPatchVectorField
        (
            const fixedPathDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedPathDisplacementFvPatchVectorField
        (
            const fixedPathDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedPathDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedPathDisplacementFvPatchVectorField
        (
            const fixedPathDisplacementFvPatchVectorField&,
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
                new fixedPathDisplacementFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Access functions

        // //- Update the coefficients associated with the patch field
        //     virtual vectorField& totalDisp()
        //     {
        //         return totalDisp_;
        //     }

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
