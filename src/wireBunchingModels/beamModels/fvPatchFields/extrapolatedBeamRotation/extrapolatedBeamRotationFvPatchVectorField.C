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

\*---------------------------------------------------------------------------*/

#include "extrapolatedBeamRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

extrapolatedBeamRotationFvPatchVectorField::
extrapolatedBeamRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    timeIndex_(-1)
{
    fvPatchVectorField::operator==(patchInternalField());
}


extrapolatedBeamRotationFvPatchVectorField::
extrapolatedBeamRotationFvPatchVectorField
(
    const extrapolatedBeamRotationFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(tdpvf, p, iF, mapper),
    timeIndex_(-1)
{}


extrapolatedBeamRotationFvPatchVectorField::
extrapolatedBeamRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    timeIndex_(-1)
{
    // if (dict.found("value"))
    // {
    //     Field<vector>::operator==(vectorField("value", dict, p.size()));
    // }
    // else
    // {
    //     fvPatchVectorField::operator==(patchInternalField());
    // }
}


extrapolatedBeamRotationFvPatchVectorField::
extrapolatedBeamRotationFvPatchVectorField
(
    const extrapolatedBeamRotationFvPatchVectorField& tdpvf
)
:
    fixedValueFvPatchVectorField(tdpvf),
    timeIndex_(-1)
{}


extrapolatedBeamRotationFvPatchVectorField::
extrapolatedBeamRotationFvPatchVectorField
(
    const extrapolatedBeamRotationFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tdpvf, iF),
    timeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void extrapolatedBeamRotationFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void extrapolatedBeamRotationFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    // const extrapolatedBeamRotationFvPatchVectorField& dmptf =
    //     refCast<const extrapolatedBeamRotationFvPatchVectorField>(ptf);
}


// Update the coefficients associated with the patch field
void extrapolatedBeamRotationFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (timeIndex_ < this->db().time().timeIndex())
    {
        timeIndex_ = this->db().time().timeIndex();

        vectorField theta (this->patchInternalField());
        forAll(theta, faceI)
        {
            // Prevent torsion
            theta[faceI].x() = 0;

            // Prevent bending
            // theta[faceI].y() = 0;
            // theta[faceI].z() = 0;
        }

        // Info << theta << endl;

        fvPatchField<vector>::operator==(theta);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

// Write
void extrapolatedBeamRotationFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, extrapolatedBeamRotationFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
