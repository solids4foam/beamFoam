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

#include "fixedPathRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
// #include "fvcMeshPhi.H"
// #include "pointMesh.H"
// #include "pointFields.H"
// #include "fixedValuePointPatchFields.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedPathRotationFvPatchVectorField::
fixedPathRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    fileName_(),
    pathPtr_(),
    velocity_(0)
{}


fixedPathRotationFvPatchVectorField::
fixedPathRotationFvPatchVectorField
(
    const fixedPathRotationFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    fileName_(),
    pathPtr_
    (
        new HermiteSpline
        (
            ptf.pathPtr_().points(),
            ptf.pathPtr_().tangents()
        )
    ),
    velocity_(0)
{}


fixedPathRotationFvPatchVectorField::
fixedPathRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    fileName_(dict.lookup("pathFileName")),
    pathPtr_(),
    velocity_(readScalar(dict.lookup("velocity")))
{
    IFstream ifs(fileName_);

    vectorField points(ifs);    
    vectorField tangents(ifs);    
    
    // pathPtr_.set(new HermiteSpline(ifs));
    pathPtr_.set(new HermiteSpline(points, tangents));
    
    // Info << pathPtr_().points() << endl;
    // Info << pathPtr_().tangents() << endl;
}


fixedPathRotationFvPatchVectorField::
fixedPathRotationFvPatchVectorField
(
    const fixedPathRotationFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    fileName_(pivpvf.fileName_),
    pathPtr_
    (
        new HermiteSpline
        (
            pivpvf.pathPtr_().points(),
            pivpvf.pathPtr_().tangents()
        )
    ),
    velocity_(pivpvf.velocity_)
{}


fixedPathRotationFvPatchVectorField::
fixedPathRotationFvPatchVectorField
(
    const fixedPathRotationFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    fileName_(pivpvf.fileName_),
    pathPtr_
    (
        new HermiteSpline
        (
            pivpvf.pathPtr_().points(),
            pivpvf.pathPtr_().tangents()
        )
    ),
    velocity_(pivpvf.velocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedPathRotationFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedPathRotationFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    // const fixedPathRotationFvPatchVectorField& dmptf =
    //     refCast<const fixedPathRotationFvPatchVectorField>(ptf);
    // totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void fixedPathRotationFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    scalar curArcLength = velocity_*this->db().time().timeOutputValue();

    vectorField tangent
    (
        patch().size(),
        pathPtr_().firstDerivative(curArcLength)
    );
    tangent /= mag(tangent) + SMALL;

    vector i(1, 0, 0);
    
    tensor R = rotationTensor(i, tangent[0]);

    vector rotAxis
    (
        R.zy() - R.yz(),
        R.xz() - R.zx(),
        R.yx() - R.xy()
    );
    
    scalar magRotAngle = asin(mag(rotAxis)/2);
    
    rotAxis /= mag(rotAxis);

    vector rotAngle = magRotAngle*rotAxis;

    vectorField theta(patch().size(), rotAngle);
    // Info << rotAngle*180/M_PI << endl;

    // const fvMesh& mesh = patch().boundaryMesh().mesh();

    // if
    // (
    //     mesh.moving() // Updated lagrangian
    // )
    // {
    //     const volVectorField& totW =
    //         db().lookupObject<volVectorField>("totW");

    //     disp -= totW.boundaryField()[patch().index()];
    // }
    
    fvPatchField<vector>::operator==(theta);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void fixedPathRotationFvPatchVectorField::write(Ostream& os) const
{
    // if (dispSeries_.size())
    // {
    //     os.writeKeyword("displacementSeries") << nl;
    //     os << token::BEGIN_BLOCK << nl;
    //     dispSeries_.write(os);
    //     os << token::END_BLOCK << nl;
    // }

    os.writeKeyword("pathFileName")
        << fileName_ << token::END_STATEMENT << nl;
    
    os.writeKeyword("velocity")
        << velocity_ << token::END_STATEMENT << nl;

    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedPathRotationFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
