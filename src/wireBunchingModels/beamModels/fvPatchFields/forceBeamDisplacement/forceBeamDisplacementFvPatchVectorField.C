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

#include "forceBeamDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

forceBeamDisplacementFvPatchVectorField::
forceBeamDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    forceSeries_()
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


forceBeamDisplacementFvPatchVectorField::
forceBeamDisplacementFvPatchVectorField
(
    const forceBeamDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    force_(tdpvf.force_, mapper),
    forceSeries_(tdpvf.forceSeries_)
{}


forceBeamDisplacementFvPatchVectorField::
forceBeamDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    // force_("force", dict, p.size()),
    forceSeries_()
{
    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }
  
    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }
    
    // fvPatchVectorField::operator=(patchInternalField());
    // gradient() = vector::zero;
    
    // Check if traction is time-varying
    if (dict.found("forceSeries"))
    {
        Info<< "force is time-varying" << endl;
        forceSeries_ =
            interpolationTable<vector>(dict.subDict("forceSeries"));
    }
    else
    {
        force_ = vectorField("force", dict, p.size());
    }
}


forceBeamDisplacementFvPatchVectorField::
forceBeamDisplacementFvPatchVectorField
(
    const forceBeamDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    force_(tdpvf.force_),
    forceSeries_(tdpvf.forceSeries_)
{}


forceBeamDisplacementFvPatchVectorField::
forceBeamDisplacementFvPatchVectorField
(
    const forceBeamDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    force_(tdpvf.force_),
    forceSeries_(tdpvf.forceSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void forceBeamDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    force_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void forceBeamDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);
    
    const forceBeamDisplacementFvPatchVectorField& dmptf =
        refCast<const forceBeamDisplacementFvPatchVectorField>(ptf);
    
    force_.rmap(dmptf.force_, addr);
}


// Update the coefficients associated with the patch field
void forceBeamDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Info << "forceBeamDisplacementFvPatchVectorField::updateCoeffs()"
    //      << endl;
    
    if (forceSeries_.size())
    {
        force_ = forceSeries_(this->db().time().timeOutputValue());
    }

    if (dimensionedInternalField().name() == "DW")
    {
        const vectorField DTheta =
            patch().lookupPatchField<volVectorField, vector>("DTheta");

        const tensorField CQDTheta =
            patch().lookupPatchField<surfaceTensorField, tensor>("CQDTheta");

        const tensorField CQDW =
            patch().lookupPatchField<surfaceTensorField, tensor>("CQDW");

        const vectorField explicitQ =
            patch().lookupPatchField<surfaceVectorField, vector>("explicitQ");
        // const vectorField explicitQ =
        //     patch().lookupPatchField<surfaceVectorField, vector>("Q_0");

        forAll(force_, faceI)
        {
            tensor invCQDW = inv(CQDW[faceI]);

            gradient()[faceI] =
            (
                invCQDW
              & (
                    force_[faceI] - explicitQ[faceI]
                  - (CQDTheta[faceI] & DTheta[faceI])
                )
            );
        }
    }
    else
    {
        const vectorField Theta =
            patch().lookupPatchField<volVectorField, vector>("Theta");

        const tensorField CQTheta =
            patch().lookupPatchField<surfaceTensorField, tensor>("CQTheta");

        const tensorField CQW =
            patch().lookupPatchField<surfaceTensorField, tensor>("CQW");

        const vectorField explicitQ =
            patch().lookupPatchField<surfaceVectorField, vector>("explicitQ");

        forAll(force_, faceI)
        {
            tensor invCQW = inv(CQW[faceI]);

            gradient()[faceI] =
            (
                invCQW
              & (
                    force_[faceI] - explicitQ[faceI]
                  - (CQTheta[faceI] & Theta[faceI])
                )
            );
        }
    }
    
    // fixedGradientFvPatchVectorField::updateCoeffs();
}


// Write
void forceBeamDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    
    if (forceSeries_.size())
    {
        os.writeKeyword("forceSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        forceSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    // else
    // {
        force_.writeEntry("force", os);
    // }
    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, forceBeamDisplacementFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
