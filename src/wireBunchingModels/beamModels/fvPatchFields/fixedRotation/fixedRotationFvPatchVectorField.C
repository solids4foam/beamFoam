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

#include "fixedRotationFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedRotationFvPatchVectorField::
fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    totalTheta_(p.size(), vector::zero),
    thetaSeries_(),
    timeIndex_(-1)
{
   // fvPatchVectorField::operator==(patchInternalField());
}


fixedRotationFvPatchVectorField::
fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(tdpvf, p, iF, mapper),
    totalTheta_(tdpvf.totalTheta_, mapper),
    thetaSeries_(tdpvf.thetaSeries_),
    timeIndex_(-1)
{}


fixedRotationFvPatchVectorField::
fixedRotationFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    totalTheta_("value", dict, p.size()),
    thetaSeries_(),
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
    
    // Check if rotation is time-varying
    if (dict.found("thetaSeries"))
    {
        Info<< "   Rotation is time-varying" << endl;
        thetaSeries_ =
            interpolationTable<vector>(dict.subDict("thetaSeries"));

        fvPatchField<vector>::operator==
        (
            thetaSeries_(this->db().time().timeOutputValue())
        );
    }
}


fixedRotationFvPatchVectorField::
fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& tdpvf
)
:
    fixedValueFvPatchVectorField(tdpvf),
    totalTheta_(tdpvf.totalTheta_),
    thetaSeries_(tdpvf.thetaSeries_),
    timeIndex_(-1)
{}


fixedRotationFvPatchVectorField::
fixedRotationFvPatchVectorField
(
    const fixedRotationFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(tdpvf, iF),
    totalTheta_(tdpvf.totalTheta_),
    thetaSeries_(tdpvf.thetaSeries_),
    timeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedRotationFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    
    totalTheta_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedRotationFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
    
    // SB uncommented this line
    const fixedRotationFvPatchVectorField& dmptf =
         refCast<const fixedRotationFvPatchVectorField>(ptf);
         
     totalTheta_.rmap(dmptf.totalTheta_, addr);
}


// Update the coefficients associated with the patch field
void fixedRotationFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    vectorField theta = totalTheta_;

    if (thetaSeries_.size())
    {
        theta = thetaSeries_(this->db().time().timeOutputValue());
    }
    
    // SB: Dont know how theta will be updated for Updated Lagrangian Formulation.
    // Check with ZT
    
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    
    if(mesh.moving())
    {
    	//const volVectorField& DTheta =
           // db().lookupObject<volVectorField>("DTheta");

      // theta = theta - thetaCurr.boundaryField()[patch().index()];
      
     // theta += DTheta;
        
      //  theta -= thetaCurr;
        
    }

/*
    if (timeIndex_ < this->db().time().timeIndex())
    {
        timeIndex_ = this->db().time().timeIndex();
        
        vectorField theta = this->patchInternalField();
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
    
*/
    fvPatchField<vector>::operator==(theta);
    fixedValueFvPatchVectorField::updateCoeffs();
}

// Write
void fixedRotationFvPatchVectorField::write(Ostream& os) const
{
	if (thetaSeries_.size())
    {
        os.writeKeyword("thetaSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        thetaSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    
    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, fixedRotationFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
