/*---------------------------------------------------------------------------* \
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

#include "axialForceTransverseDisplacementNRFvPatchVectorField.H"
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

axialForceTransverseDisplacementNRFvPatchVectorField::
axialForceTransverseDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refDisp_(p.size(), vector::zero),
    refDispSeries_(),
    axialForce_(p.size(), 0),
    axialForceSeries_()
{}


axialForceTransverseDisplacementNRFvPatchVectorField::
axialForceTransverseDisplacementNRFvPatchVectorField
(
    const axialForceTransverseDisplacementNRFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    refDisp_(ptf.refDisp_, mapper),
    refDispSeries_(ptf.refDispSeries_),
    axialForce_(ptf.axialForce_, mapper),
    axialForceSeries_(ptf.axialForceSeries_)
{}


axialForceTransverseDisplacementNRFvPatchVectorField::
axialForceTransverseDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    refDisp_(p.size(), vector::zero),
    // refDisp_("refDisp", dict, p.size()),
    refDispSeries_(),
    axialForce_(p.size(), 0),
    // axialForce_("axialForce", dict, p.size()),
    axialForceSeries_()
{
    if (dict.found("value"))
    {
        fixedValueFvPatchVectorField::operator==
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchVectorField::operator==
        (
            patchInternalField()
        );
    }
    
    // Check if axial force is time-varying
    if (dict.found("axialForceSeries"))
    {
        Info<< "    axial force is time-varying" << endl;
        axialForceSeries_ =
            interpolationTable<scalar>(dict.subDict("axialForceSeries"));

        axialForce_ = axialForceSeries_(this->db().time().timeOutputValue());
    }
    else
    {
        axialForce_ = scalarField("axialForce", dict, p.size());
    }
    
    // Check if refDisp is time-varying
    if (dict.found("refDispSeries"))
    {
        Info<< "    refDisp is time-varying" << endl;
        refDispSeries_ =
            interpolationTable<vector>(dict.subDict("refDispSeries"));

        refDisp_ = refDispSeries_(this->db().time().timeOutputValue());
    }
    else
    {
        refDisp_ = vectorField("refDisp", dict, p.size());
    }
}


axialForceTransverseDisplacementNRFvPatchVectorField::
axialForceTransverseDisplacementNRFvPatchVectorField
(
    const axialForceTransverseDisplacementNRFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    refDisp_(pivpvf.refDisp_),
    refDispSeries_(pivpvf.refDispSeries_),
    axialForce_(pivpvf.axialForce_),
    axialForceSeries_(pivpvf.axialForceSeries_)
{}


axialForceTransverseDisplacementNRFvPatchVectorField::
axialForceTransverseDisplacementNRFvPatchVectorField
(
    const axialForceTransverseDisplacementNRFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    refDisp_(pivpvf.refDisp_),
    refDispSeries_(pivpvf.refDispSeries_),
    axialForce_(pivpvf.axialForce_),
    axialForceSeries_(pivpvf.axialForceSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void axialForceTransverseDisplacementNRFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    refDisp_.autoMap(m);
    axialForce_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void axialForceTransverseDisplacementNRFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const axialForceTransverseDisplacementNRFvPatchVectorField& dmptf =
        refCast<const axialForceTransverseDisplacementNRFvPatchVectorField>(ptf);

    refDisp_.rmap(dmptf.refDisp_, addr);
    axialForce_.rmap(dmptf.axialForce_, addr);
}


void axialForceTransverseDisplacementNRFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (axialForceSeries_.size())
    {
        axialForce_ = axialForceSeries_(this->db().time().timeOutputValue());
    }

    if (refDispSeries_.size())
    {
        refDisp_ = refDispSeries_(this->db().time().timeOutputValue());
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    
    if
    (
        mesh.moving() // Updated lagrangian
    )
    {
        if (refDispSeries_.size())
        {
            refDisp_ -=
                refDispSeries_
                (
                    this->db().time().timeOutputValue()
                  - this->db().time().deltaT().value()
                );

            Info << refDisp_ << endl;
        }
        else
        {        
            const volVectorField& totW =
                db().lookupObject<volVectorField>("totW");
            refDisp_ -= totW.boundaryField()[patch().index()];
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void axialForceTransverseDisplacementNRFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchField<vector>& DW =
        const_cast<fvPatchField<vector>& >
        (
            patch().lookupPatchField<volVectorField, vector>("DW")
        );

    const fvPatchField<vector>& DTheta =
        patch().lookupPatchField<volVectorField, vector>("DTheta");

    const tensorField& Lambda =
        patch().lookupPatchField<surfaceTensorField, tensor>("Lambda");

    const vectorField& dR0Ds =
        patch().lookupPatchField<surfaceVectorField, vector>("dR0Ds");

    vectorField t(Lambda & dR0Ds);

    const tensorField& CQW =
        patch().lookupPatchField<surfaceTensorField, tensor>("CQW");

    const tensorField& CQTheta =
        patch().lookupPatchField<surfaceTensorField, tensor>("CQTheta");

    const vectorField& explicitQ =
        patch().lookupPatchField<surfaceVectorField, vector>("explicitQ");

    const scalarField delta (1.0/patch().deltaCoeffs());
 
    tensorField invCQW (inv(CQW));

    vectorField DDW
	(
	    (
	        invCQW
	      & (
	            axialForce_*t
	          - ((t*t) & explicitQ)
	          - (((t*t) & CQTheta) & DTheta)
	        )
	    )*delta
	);	
    DDW +=
    (
        (I-(t*t))
      & (
            (refDisp_ - *this)
          - DW.patchInternalField()
        )
    );

    vectorField newDW
	(
		DW.patchInternalField() + DDW
	);
    // Info << "newDW = " << newDW << endl;

    DW = newDW;

    fixedValueFvPatchField<vector>::operator==((*this) + newDW);
    
    fixedValueFvPatchVectorField::evaluate();
}

  
void axialForceTransverseDisplacementNRFvPatchVectorField::write
(
    Ostream& os
) const
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
    
    if (refDispSeries_.size())
    {
        os.writeKeyword("refDispSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        refDispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    axialForceTransverseDisplacementNRFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
