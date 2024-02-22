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

\*---------------------------------------------------------------------------*/

#include "axialForceTransverseDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
// #include "fvcMeshPhi.H"
// #include "pointMesh.H"
// #include "pointFields.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

axialForceTransverseDisplacementFvPatchVectorField::
axialForceTransverseDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refDisp_(p.size(), vector::zero),
    axialForce_(p.size(), 0),
    axialForceSeries_()
{}


axialForceTransverseDisplacementFvPatchVectorField::
axialForceTransverseDisplacementFvPatchVectorField
(
    const axialForceTransverseDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    refDisp_(ptf.refDisp_, mapper),
    axialForce_(ptf.axialForce_, mapper),
    axialForceSeries_(ptf.axialForceSeries_)
{}


axialForceTransverseDisplacementFvPatchVectorField::
axialForceTransverseDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    refDisp_("refDisp", dict, p.size()),
    axialForce_("axialForce", dict, p.size()),
    axialForceSeries_()
{
    // Check if axial force is time-varying
    if (dict.found("axialForceSeries"))
    {
        Info<< "    axial force is time-varying" << endl;
        axialForceSeries_ =
            interpolationTable<scalar>(dict.subDict("axialForceSeries"));

        axialForce_ = axialForceSeries_(this->db().time().timeOutputValue());
    }
}


axialForceTransverseDisplacementFvPatchVectorField::
axialForceTransverseDisplacementFvPatchVectorField
(
    const axialForceTransverseDisplacementFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    refDisp_(pivpvf.refDisp_),
    axialForce_(pivpvf.axialForce_),
    axialForceSeries_(pivpvf.axialForceSeries_)
{}


axialForceTransverseDisplacementFvPatchVectorField::
axialForceTransverseDisplacementFvPatchVectorField
(
    const axialForceTransverseDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    refDisp_(pivpvf.refDisp_),
    axialForce_(pivpvf.axialForce_),
    axialForceSeries_(pivpvf.axialForceSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void axialForceTransverseDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    refDisp_.autoMap(m);
    axialForce_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void axialForceTransverseDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const axialForceTransverseDisplacementFvPatchVectorField& dmptf =
        refCast<const axialForceTransverseDisplacementFvPatchVectorField>(ptf);

    refDisp_.rmap(dmptf.refDisp_, addr);
    axialForce_.rmap(dmptf.axialForce_, addr);
}


void axialForceTransverseDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // vectorField disp = refDisp_;

    if (axialForceSeries_.size())
    {
        axialForce_ = axialForceSeries_(this->db().time().timeOutputValue());
    }

    // if (dispSeries_.size())
    // {
    //     disp = dispSeries_(this->db().time().timeOutputValue());
    // }

    // if (dimensionedInternalField().name() == "DW")
    // {
    //     // Incremental approach, so we will set the increment of displacement
    //     // Lookup the old displacement field and subtract it from the ref
    //     // displacement
    //     const volVectorField& Wold =
    //         db().lookupObject<volVectorField>("W").oldTime();

    //     disp -= Wold.boundaryField()[patch().index()];
    // }


    if (internalField().name() == "DW")
    {
        const tensorField DLambda =
            patch().lookupPatchField<surfaceTensorField, tensor>("DLambda");

        const vectorField dRuDs =
            patch().lookupPatchField<surfaceVectorField, tensor>("dRuDs");

        const tensorField CQDTheta =
            patch().lookupPatchField<surfaceTensorField, vector>("CQDTheta");

        const tensorField CQDW =
            patch().lookupPatchField<surfaceTensorField, vector>("CQDW");

        const vectorField oldForce =
            patch().lookupPatchField<surfaceVectorField, vector>("Q_0");

        const vectorField dQ =
            patch().lookupPatchField<surfaceVectorField, vector>("dQ");

        const vectorField DTheta =
            patch().lookupPatchField<volVectorField, vector>("DTheta");

        const scalarField deltaCoeffs = patch().deltaCoeffs();

        vectorField t ((DLambda & dRuDs));

        // Current transverse force
        vectorField Qt
	(
            (CQDW & (refDisp_ - patchInternalField()))*deltaCoeffs
          + (CQDTheta & DTheta) + oldForce + dQ
	);
	  Qt += t*axialForce_ - ((t*t) & Qt);

        // Current displacement
        tensorField invCQDW (inv(CQDW));
        vectorField disp
	(
            patchInternalField()
          + (
                invCQDW
              & (
                    Qt
                  - (CQDTheta & DTheta)
                  - oldForce - dQ
                )
            )/deltaCoeffs
    	);
        fvPatchField<vector>::operator==(disp);
    }
    else
    {
        const tensorField Lambda =
            patch().lookupPatchField<surfaceTensorField, tensor>("Lambda");

        const vectorField dR0Ds =
            patch().lookupPatchField<surfaceVectorField, vector>("dR0Ds");

        // const vectorField i =
        //     patch().lookupPatchField<surfaceVectorField, vector>("i");

        const tensorField CQTheta =
            patch().lookupPatchField<surfaceTensorField, tensor>("CQTheta");

        const tensorField CQW =
            patch().lookupPatchField<surfaceTensorField, tensor>("CQW");

        const vectorField explicitQ =
            patch().lookupPatchField<surfaceVectorField, vector>("explicitQ");

        const vectorField Theta =
            patch().lookupPatchField<volVectorField, vector>("Theta");

        const scalarField deltaCoeffs = patch().deltaCoeffs();

        vectorField t ((Lambda & dR0Ds));

        // Current transverse force
        vectorField Qt
	(
            (CQW & (refDisp_ - patchInternalField()))*deltaCoeffs
          + (CQTheta & Theta) + explicitQ
	);
        Qt += t*axialForce_ - ((t*t) & Qt);

        // Current displacement
        tensorField invCQW (inv(CQW));
        vectorField disp
	(
            patchInternalField()
          + (
                invCQW
              & (
                    Qt - (CQTheta & Theta) - explicitQ
                )
            )/deltaCoeffs
	);
        fvPatchField<vector>::operator==(disp);
    }

    // fixedValueFvPatchVectorField::updateCoeffs();
}

void axialForceTransverseDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    refDisp_.writeEntry("refDisp", os);
    axialForce_.writeEntry("axialForce", os);
    if (axialForceSeries_.size())
    {
        os.writeKeyword("axialForceSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        axialForceSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    axialForceTransverseDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
