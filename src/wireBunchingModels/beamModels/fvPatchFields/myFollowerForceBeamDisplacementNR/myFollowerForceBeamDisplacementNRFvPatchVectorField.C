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

#include "myFollowerForceBeamDisplacementNRFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myFollowerForceBeamDisplacementNRFvPatchVectorField::
myFollowerForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    followerPlaneNormal_(p.size(), vector::zero),
    followerForce_(p.size(), 0),
    followerForceDir_(p.size(), vector::zero),
    outOfPlaneForce_(p.size(), 0),
    followerForceSeries_(),
    outOfPlaneForceSeries_(),
    timeIndex_(-1)
{
    // fixedFvPatchVectorField::operator=(patchInternalField());
    // gradient() = vector::zero;
}


myFollowerForceBeamDisplacementNRFvPatchVectorField::
myFollowerForceBeamDisplacementNRFvPatchVectorField
(
    const myFollowerForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, p, iF, mapper),
    followerPlaneNormal_(tdpvf.followerPlaneNormal_, mapper),
    followerForce_(tdpvf.followerForce_, mapper),
    followerForceDir_(tdpvf.followerForceDir_, mapper),
    outOfPlaneForce_(tdpvf.outOfPlaneForce_, mapper),
    followerForceSeries_(tdpvf.followerForceSeries_),
    outOfPlaneForceSeries_(tdpvf.outOfPlaneForceSeries_),
    timeIndex_(tdpvf.timeIndex_)
{}


myFollowerForceBeamDisplacementNRFvPatchVectorField::
myFollowerForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    // forceBeamDisplacementNRFvPatchVectorField(p, iF, dict),
    followerPlaneNormal_("followerPlaneNormal", dict, p.size()),
    followerForce_(p.size(), 0),
    followerForceDir_(p.size(), vector::zero),
    outOfPlaneForce_(p.size(), 0),
    followerForceSeries_(),
    outOfPlaneForceSeries_(),
    timeIndex_(-1)
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

    if (dict.found("myFollowerForceSeries"))
    {
        Info<< "follower force is time-varying" << endl;
        followerForceSeries_ =
            interpolationTable<scalar>(dict.subDict("myFollowerForceSeries"));

        followerForce_ =
            followerForceSeries_(this->db().time().timeOutputValue());
    }
    else
    {
        FatalErrorIn
        (
            "myFollowerForceBeamDisplacementNRFvPatchVectorField::"
            "myFollowerForceBeamDisplacementNRFvPatchVectorField"
            "("
            "    const fvPatch& p,"
            "    const DimensionedField<vector, volMesh>& iF,"
            "    const dictionary& dict"
            ")"
        )
            << "myFollowerForceSeries is not defined"
            << abort(FatalError);
    }

    if (dict.found("outOfPlaneForceSeries"))
    {
        Info<< "out-to-plane force is time-varying" << endl;
        outOfPlaneForceSeries_ =
            interpolationTable<scalar>(dict.subDict("outOfPlaneForceSeries"));

        outOfPlaneForce_ =
            outOfPlaneForceSeries_(this->db().time().timeOutputValue());
    }
    // else
    // {
    //     FatalErrorIn
    //     (
    //         "followerForceBeamDisplacementNRFvPatchVectorField::"
    //         "followerForceBeamDisplacementNRFvPatchVectorField"
    //         "("
    //         "    const fvPatch& p,"
    //         "    const DimensionedField<vector, volMesh>& iF,"
    //         "    const dictionary& dict"
    //         ")"
    //     )
    //         << "outOfPlaneForceSeries is not defined"
    //         << abort(FatalError);
    // }
}


myFollowerForceBeamDisplacementNRFvPatchVectorField::
myFollowerForceBeamDisplacementNRFvPatchVectorField
(
    const myFollowerForceBeamDisplacementNRFvPatchVectorField& tdpvf
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf),
    followerPlaneNormal_(tdpvf.followerPlaneNormal_),
    followerForce_(tdpvf.followerForce_),
    followerForceDir_(tdpvf.followerForceDir_),
    outOfPlaneForce_(tdpvf.outOfPlaneForce_),
    followerForceSeries_(tdpvf.followerForceSeries_),
    outOfPlaneForceSeries_(tdpvf.outOfPlaneForceSeries_),
    timeIndex_(tdpvf.timeIndex_)
{}


myFollowerForceBeamDisplacementNRFvPatchVectorField::
myFollowerForceBeamDisplacementNRFvPatchVectorField
(
    const myFollowerForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, iF),
    followerPlaneNormal_(tdpvf.followerPlaneNormal_),
    followerForce_(tdpvf.followerForce_),
    followerForceDir_(tdpvf.followerForceDir_),
    outOfPlaneForce_(tdpvf.outOfPlaneForce_),
    followerForceSeries_(tdpvf.followerForceSeries_),
    outOfPlaneForceSeries_(tdpvf.outOfPlaneForceSeries_),
    timeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myFollowerForceBeamDisplacementNRFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    forceBeamDisplacementNRFvPatchVectorField::autoMap(m);
    followerPlaneNormal_.autoMap(m);
    followerForce_.autoMap(m);
    outOfPlaneForce_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void myFollowerForceBeamDisplacementNRFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    forceBeamDisplacementNRFvPatchVectorField::rmap(ptf, addr);
    
    const myFollowerForceBeamDisplacementNRFvPatchVectorField& dmptf =
        refCast<const myFollowerForceBeamDisplacementNRFvPatchVectorField>(ptf);

    followerPlaneNormal_.rmap(dmptf.followerPlaneNormal_, addr);
    followerForce_.rmap(dmptf.followerForce_, addr);
    outOfPlaneForce_.rmap(dmptf.outOfPlaneForce_, addr);
    timeIndex_ = dmptf.timeIndex_;
}


// Update the coefficients associated with the patch field
void myFollowerForceBeamDisplacementNRFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if (timeIndex_ < this->db().time().timeIndex())
    {
        timeIndex_ = this->db().time().timeIndex();

        followerForce_ = followerForceSeries_
        (
            this->db().time().timeOutputValue()
        );

        if (outOfPlaneForceSeries_.size())
        {
            outOfPlaneForce_ = outOfPlaneForceSeries_
            (
                this->db().time().timeOutputValue()
            );
        }

//        const vectorField& dR0Ds =
//            patch().lookupPatchField<surfaceVectorField, vector>("dR0Ds");

//        const vectorField dWdS =
//            db().lookupObject<volVectorField>("W")
//           .boundaryField()[this->patch().index()].snGrad();

        // const vectorField dWdS =
        //     db().lookupObject<volVectorField>("W").oldTime()
        //    .boundaryField()[this->patch().index()].snGrad();

//        vectorField dRdS = dR0Ds + dWdS;
//        vectorField t = dRdS/(mag(dRdS) + SMALL);

	    const tensorField& Lambda =
            patch().lookupPatchField<surfaceTensorField, tensor>("Lambda");
/*
        if ( timeIndex_ < 51)
        {
            followerForceDir_ = (followerPlaneNormal_ ^ t);
        }
*/
//		followerForceDir_ = (followerPlaneNormal_ ^ t);


		const scalarField tx = Lambda.component(tensor::XZ);
		const scalarField ty = Lambda.component(tensor::YZ);
		const scalarField tz = Lambda.component(tensor::ZZ);
		
		const vector i(1, 0, 0);
		const vector j(0, 1, 0);
		const vector k(0, 0, 1);
		
		const vectorField tnew = (i*tx + j*ty + k*tz);
		Info << "tnew " << tnew << endl;
		
		const vectorField t = tnew/mag(tnew);
		
//		tnew = tnew/magt;
		
//		Info << "tnew normalised: " << tnew << endl;
//		const vectorField tnew = (tx, ty, tz);
//		t = t/mag(t);
		
		followerForceDir_ = t;
      
        Info << "Lambda " << Lambda << endl;
//        Info << "tnew " << tnew << endl;
        Info << "followerForceDir: " << followerForceDir_ << endl;
        Info << "followerForce: " << followerForce_ << endl;
        
        vectorField curForce =
            followerForce_*followerForceDir_;
          // + outOfPlaneForce_*followerPlaneNormal_;
        
        Info << "curForce: " << curForce << endl;
        force() = curForce;
    }

    // fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void myFollowerForceBeamDisplacementNRFvPatchVectorField::write(Ostream& os) const
{
    forceBeamDisplacementNRFvPatchVectorField::write(os);
    
    if (followerForceSeries_.size())
    {
        os.writeKeyword("followerForceSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        followerForceSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    if (outOfPlaneForceSeries_.size())
    {
        os.writeKeyword("outOfPlaneForceSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        outOfPlaneForceSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    followerPlaneNormal_.writeEntry("followerPlaneNormal", os);
    followerForce_.writeEntry("followerForce", os);
    outOfPlaneForce_.writeEntry("outOfPlaneForce", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    myFollowerForceBeamDisplacementNRFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
