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
    Foam::axialForceTransverseDisplacementNRFvPatchVectorField

Description
    Fixed transverse displacement and axial force boundary condition for beams

SourceFiles
    axialForceTransverseDisplacementNRFvPatchVectorField.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef axialForceTransverseDisplacementNRFvPatchVectorField_H
#define axialForceTransverseDisplacementNRFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class axialForceTransverseDisplacementNRFvPatch Declaration
\*---------------------------------------------------------------------------*/

class axialForceTransverseDisplacementNRFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Reference displacement
        vectorField refDisp_;

        //- Displacement time series
        interpolationTable<vector> refDispSeries_;

        //- Axial force
        scalarField axialForce_;

        //- Displacement time series
        interpolationTable<scalar> axialForceSeries_;

    // Private Functions

public:

    //- Runtime type information
    TypeName("axialForceTransverseDisplacementNR");


    // Constructors

        //- Construct from patch and internal field
        axialForceTransverseDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        axialForceTransverseDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  axialForceTransverseDisplacementNRFvPatchVectorField
        //  onto a new patch
        axialForceTransverseDisplacementNRFvPatchVectorField
        (
            const axialForceTransverseDisplacementNRFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        axialForceTransverseDisplacementNRFvPatchVectorField
        (
            const axialForceTransverseDisplacementNRFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new axialForceTransverseDisplacementNRFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        axialForceTransverseDisplacementNRFvPatchVectorField
        (
            const axialForceTransverseDisplacementNRFvPatchVectorField&,
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
                new axialForceTransverseDisplacementNRFvPatchVectorField
                (*this, iF)
            );
        }



    // Member functions

        // Access functions

            //- Return reference displacement
            const vectorField& refDisp() const
            {
                return refDisp_;
            }

            //- Return reference displacement
            vectorField& refDisp()
            {
                return refDisp_;
            }

            //- Return axial force
            const scalarField& axialForce() const
            {
                return axialForce_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
            );

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
