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

#include "followerDisplacementFvPatchVectorField.H"
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

followerDisplacementFvPatchVectorField::followerDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    followerDispSeries_()
{}


followerDisplacementFvPatchVectorField::followerDisplacementFvPatchVectorField
(
    const followerDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    totalDisp_(ptf.totalDisp_, mapper),
    followerDispSeries_(ptf.followerDispSeries_)
{}


followerDisplacementFvPatchVectorField::followerDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    totalDisp_("value", dict, p.size()),
    followerDispSeries_()
{
    // Check if displacement is time-varying
    if (dict.found("followerDisplacementSeries"))
    {
        Info<< " follower displacement is time-varying" << endl;
        followerDispSeries_ =
            interpolationTable<vector>(dict.subDict("followerDisplacementSeries"));

        fvPatchField<vector>::operator==
        (
            followerDispSeries_(this->db().time().timeOutputValue())
        );
    }
}


followerDisplacementFvPatchVectorField::followerDisplacementFvPatchVectorField
(
    const followerDisplacementFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    totalDisp_(pivpvf.totalDisp_),
    followerDispSeries_(pivpvf.followerDispSeries_)
{}


followerDisplacementFvPatchVectorField::followerDisplacementFvPatchVectorField
(
    const followerDisplacementFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    totalDisp_(pivpvf.totalDisp_),
    followerDispSeries_(pivpvf.followerDispSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void followerDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

    totalDisp_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void followerDisplacementFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const followerDisplacementFvPatchVectorField& dmptf =
        refCast<const followerDisplacementFvPatchVectorField>(ptf);

    totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void followerDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField disp = totalDisp_;

    if (followerDispSeries_.size())
    {
        disp = followerDispSeries_(this->db().time().timeOutputValue());
	
	const tensorField& Lambda =
            patch().lookupPatchField<surfaceTensorField, tensor>("Lambda");
	
	// Extraction of column vectors from the rotation matrix
	// to construct the body-attached coordinate
	// frame t1,t2,t3
	
	// 1st column
	const scalarField txx = Lambda.component(tensor::XX);
	const scalarField tyx = Lambda.component(tensor::YX);
	const scalarField tzx = Lambda.component(tensor::ZX);
	
	// 2nd Column
	const scalarField txy = Lambda.component(tensor::XY);
	const scalarField tyy = Lambda.component(tensor::YY);
	const scalarField tzy = Lambda.component(tensor::ZY);
	
	// 3rd Column
	const scalarField txz = Lambda.component(tensor::XZ);
	const scalarField tyz = Lambda.component(tensor::YZ);
	const scalarField tzz = Lambda.component(tensor::ZZ);
	
	// Fixed reference coordinate system	
	const vector i(1, 0, 0);
	const vector j(0, 1, 0);
	const vector k(0, 0, 1);
	
	// Body-attached coordinate frame
	vectorField t1 = (i*txx + j*tyx + k*tzx);
	vectorField t2 = (i*txy + j*tyy + k*tzy);	
	vectorField t3 = (i*txz + j*tyz + k*tzz);
	
	// Follower displacement vector
	vectorField curFollowerDisp =   disp.component(0)*t1
					+ disp.component(1)*t2
					+ disp.component(2)*t3;
	
	disp = curFollowerDisp;
	
	Info << "Displacement " << disp << endl;
    }

    const fvMesh& mesh = patch().boundaryMesh().mesh();
    
    // if
    // (
    //     (dimensionedInternalField().name() == "DW")
    // )
    // {
    //     // Incremental approach, so we will set the increment of displacement
    //     // Lookup the old displacement field and subtract it from the total
    //     // displacement
    //     const volVectorField& Wold =
    //         db().lookupObject<volVectorField>("W").oldTime();

    //     disp -= Wold.boundaryField()[patch().index()];
    // }
    // else if
    
    if
    (
        mesh.moving() // Updated lagrangian
    )
    {
        const volVectorField& totW =
            db().lookupObject<volVectorField>("totW");

        disp -= totW.boundaryField()[patch().index()];
    }

    // Info << patch().name() << ", " << disp << endl;
    
    fvPatchField<vector>::operator==(disp);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void followerDisplacementFvPatchVectorField::write(Ostream& os) const
{
    if (followerDispSeries_.size())
    {
        os.writeKeyword("followerDisplacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        followerDispSeries_.write(os);
        os << token::END_BLOCK << nl;
    }

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    followerDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
