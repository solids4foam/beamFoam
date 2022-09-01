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

#include "momentBeamRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

momentBeamRotationFvPatchVectorField::
momentBeamRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    moment_(p.size(), vector::zero),
    momentSeries_()
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


momentBeamRotationFvPatchVectorField::
momentBeamRotationFvPatchVectorField
(
    const momentBeamRotationFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    moment_(tdpvf.moment_, mapper),
    momentSeries_(tdpvf.momentSeries_)
{}


momentBeamRotationFvPatchVectorField::
momentBeamRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    moment_(p.size(), vector::zero),
    // moment_("moment", dict, p.size()),
    momentSeries_()
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
    if (dict.found("momentSeries"))
    {
        Info<< "moment is time-varying" << endl;
        momentSeries_ =
            interpolationTable<vector>(dict.subDict("momentSeries"));
    }
    else
    {
        moment_ = vectorField("moment", dict, p.size());
    }
    
}


momentBeamRotationFvPatchVectorField::
momentBeamRotationFvPatchVectorField
(
    const momentBeamRotationFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    moment_(tdpvf.moment_),
    momentSeries_(tdpvf.momentSeries_)
{}


momentBeamRotationFvPatchVectorField::
momentBeamRotationFvPatchVectorField
(
    const momentBeamRotationFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    moment_(tdpvf.moment_),
    momentSeries_(tdpvf.momentSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void momentBeamRotationFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    moment_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void momentBeamRotationFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);
    
    const momentBeamRotationFvPatchVectorField& dmptf =
        refCast<const momentBeamRotationFvPatchVectorField>(ptf);
    
    moment_.rmap(dmptf.moment_, addr);
}


// Update the coefficients associated with the patch field
void momentBeamRotationFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Info << "void momentBeamRotationFvPatchVectorField::updateCoeffs()"
    //      << endl;
      
    if (momentSeries_.size())
    {
        moment_ = momentSeries_(this->db().time().timeOutputValue());
    }

    if (dimensionedInternalField().name() == "DTheta")
    {
        const tensorField CMDTheta =
            patch().lookupPatchField<surfaceTensorField, tensor>("CMDTheta");
 
        const tensorField CMDTheta2 =
            patch().lookupPatchField<surfaceTensorField, tensor>("CMDTheta2");

        const scalarField delta = 1.0/patch().deltaCoeffs();
    
        const vectorField explicitM =
            patch().lookupPatchField<surfaceVectorField, vector>("explicitM");        
        // const vectorField explicitM =
        //     patch().lookupPatchField<surfaceVectorField, vector>("M_0");

        const vectorField DThetaP = this->patchInternalField();

        // const tensorField invCM = inv(CMDTheta/delta);
        const tensorField invCM = inv(CMDTheta/delta + CMDTheta2);
    
        vectorField DTheta =
        (
            invCM
          & (
                moment_ - explicitM
              + (CMDTheta & DThetaP)/delta
            )
        );
    
        gradient() = (DTheta - DThetaP)/delta;
    }
    else
    {
        const tensorField CMTheta =
            patch().lookupPatchField<surfaceTensorField, tensor>("CMTheta");

        const tensorField CMTheta2 =
            patch().lookupPatchField<surfaceTensorField, tensor>("CMTheta2");
        
        const vectorField explicitM =
            patch().lookupPatchField<surfaceVectorField, vector>("explicitM");

        const scalarField delta = 1.0/patch().deltaCoeffs();
    
        const vectorField ThetaP = this->patchInternalField();

        const tensorField invCM = inv(CMTheta/delta + CMTheta2);
    
        vectorField Theta =
        (
            invCM
          & (
                moment_ - explicitM
              + (CMTheta & ThetaP)/delta
            )
        );
    
        gradient() = (Theta - ThetaP)/delta;
    }

    // fixedGradientFvPatchVectorField::updateCoeffs();
}


// Write
void momentBeamRotationFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);
    moment_.writeEntry("moment", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, momentBeamRotationFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
