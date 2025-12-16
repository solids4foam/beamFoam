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

#include "linearSpringForceBeamDisplacementNRFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupBeamModel.H"
#include "pseudoVector.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    springStiffness_(0.0),
    orthogonalStiffnessRatio_(0.0),
    springDirection_(vector::zero),
    Ks_(p.size(), tensor::zero) //,
    // extraForce_(p.size(), vector::zero),
    // forceSeries_()
{
    // fixedValueFvPatchVectorField::operator==(patchInternalField());
    // gradient() = vector::zero;
}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const linearSpringForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, p, iF, mapper),
    springStiffness_(tdpvf.springStiffness_),
    springDirection_(tdpvf.springDirection_),
    Ks_(tdpvf.Ks_, mapper) //,
    // extraForce_(tdpvf.extraForce_, mapper),
    // forceSeries_(tdpvf.forceSeries_)
{}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    forceBeamDisplacementNRFvPatchVectorField(p, iF),
    springStiffness_(0.0),
    orthogonalStiffnessRatio_(0.0),
    springDirection_(vector::zero),
    Ks_(p.size(), tensor::zero) //,
    // extraForce_(p.size(), vector::zero),
    // forceSeries_()
{
    // // NOT SURE ABOUT THIS
    // if (dict.found("value"))
    // {
    //     fixedValueFvPatchVectorField::operator==
    //     (
    //         vectorField("value", dict, p.size())
    //     );
    // }
    // else
    // {
        // fixedValueFvPatchVectorField::operator==
        // (
        //     patchInternalField()
        // );
    // }

    // Constructing the spring dictionary
    const dictionary& springDict = dict.subDict("linearSpringCoeffs");

    // Read the spring stiffness
    springStiffness_ = readScalar(springDict.lookup("springStiffness"));

    // Read the orthogonal stiffness ratio
    orthogonalStiffnessRatio_ =
          springDict.lookupOrDefault<scalar>
           (
            "orthogonalStiffnessRatio",
            0.0
           );

    
    // Read the spring direction
    springDirection_ = vector(springDict.lookup("springDirection"));

    // Normalise the spring direction to unit vector
    if (mag(springDirection_) > SMALL)
    {
        springDirection_ /= mag(springDirection_);
    }

    // Spring stiffness tensor
    Ks_ = springStiffness_*(springDirection_*springDirection_);

    // Add stiffness springs in the tangential directions to weakly enforce zero displacement in the tangential directions
    Ks_ += orthogonalStiffnessRatio_*springStiffness_*(I - sqr(springDirection_));    

    if (Pstream::master())
    {
    Info << "Spring stiffness tensor Ks = " << Ks_ << nl;
    }


    
    // if (dict.found("forceSeries"))
    // {
    //     Info<< "force is time-varying" << endl;
    //     forceSeries_ =
    //         interpolationTable<vector>(dict.subDict("forceSeries"));

    //     extraForce_ = forceSeries_(this->db().time().timeOutputValue());
    //     Info<< "Extra force" << extraForce_ << endl;
    // }

}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const linearSpringForceBeamDisplacementNRFvPatchVectorField& tdpvf
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf),
    springStiffness_(tdpvf.springStiffness_),
    springDirection_(tdpvf.springDirection_),
    Ks_(tdpvf.Ks_) //,
    // extraForce_(tdpvf.extraForce_),
    // forceSeries_(tdpvf.forceSeries_)
{}


linearSpringForceBeamDisplacementNRFvPatchVectorField::
linearSpringForceBeamDisplacementNRFvPatchVectorField
(
    const linearSpringForceBeamDisplacementNRFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    forceBeamDisplacementNRFvPatchVectorField(tdpvf, iF),
    springStiffness_(tdpvf.springStiffness_),
    springDirection_(tdpvf.springDirection_),
    Ks_(tdpvf.Ks_) //,
    // extraForce_(tdpvf.extraForce_),
    // forceSeries_(tdpvf.forceSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void linearSpringForceBeamDisplacementNRFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    forceBeamDisplacementNRFvPatchVectorField::autoMap(m);
    Ks_.autoMap(m);
    // extraForce_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void linearSpringForceBeamDisplacementNRFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    forceBeamDisplacementNRFvPatchVectorField::rmap(ptf, addr);

    const linearSpringForceBeamDisplacementNRFvPatchVectorField& dmptf =
        refCast<const linearSpringForceBeamDisplacementNRFvPatchVectorField>(ptf);
    
    Ks_.rmap(dmptf.Ks_, addr);
    // extraForce_.rmap(dmptf.extraForce_, addr);    
}


// Update the coefficients associated with the patch field
void linearSpringForceBeamDisplacementNRFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    forceBeamDisplacementNRFvPatchVectorField::updateCoeffs();
}

void linearSpringForceBeamDisplacementNRFvPatchVectorField::evaluate
(
    const Pstream::commsTypes
)
{
    // Value of displacement at the boundary
    const vectorField& pW = *this;

    // Current Spring Force calculated from previous Newton iteration
    forAll(pW, pFace)
    {
        // Spring force magnitude calculated from displacement
        // component in the spring direction
      //        scalar curSpringForceMag =
      //    springStiffness_*(pW[pFace] & springDirection_);

        // Spring force vector directly updated into base
        // forceBeamDisplacementNR class.
        // This force() is also the used in the source vector
        // of global AX = B system.

        // This force is compressive and acts opposite of spring direction
        //force()[pFace] = -curSpringForceMag*springDirection_;
          // + extraForce_[pFace];

      // ============================================================================
      // Current spring force evaluated using displacement from the previous
      // Newton–Raphson iteration.
      //
      // This force contributes directly to the RHS (B) of the global
      // linear system A * X = B solved by the beam Newton–Raphson solver.
      // ============================================================================

      // pW : patch field of beam displacement vectors (one per patch face)
      // Ks_: spring stiffness tensor (or tensor field) defining force–displacement
      //      coupling, allowing forces in arbitrary directions (not just axial)

      // ---------------------------------------------------------------------------
      // IMPORTANT:
      // The contraction (Ks_ & pW) returns a tmp<Field<vector>>.
      // We must *own* this temporary to avoid dangling references.
      // ---------------------------------------------------------------------------
      tmp<Field<vector>> tKspW = Ks_ & pW;

      // Extract the underlying Field<vector> from the tmp<>.      
      //Lifetime is guaranteed as long as tKspW remains in scope.
      const Field<vector>& KspW = tKspW();
      // Spring force is linear elastic:
      //
      //     F = -K * u
      //
      // where:
      //   K  = spring stiffness tensor
      //   u  = displacement vector
      //
      // The minus sign ensures the force is restoring (opposes displacement).
      //
      // force() is the fvPatchVectorField inherited from the base class and
      // represents the force applied by the beam restraint on each patch face.
      force()[pFace] = -KspW[pFace];
    }

    // Other necessary fields to enforce the boundary condition
    fvPatchField<vector>& DW =
        const_cast<fvPatchField<vector>&>
        (
            patch().lookupPatchField<volVectorField, vector>("DW")
        );

    const fvPatchField<vector>& DTheta =
        patch().lookupPatchField<volVectorField, vector>("DTheta");

    const tensorField& CQW =
        patch().lookupPatchField<surfaceTensorField, tensor>("CQW");

    const tensorField& CQTheta =
        patch().lookupPatchField<surfaceTensorField, tensor>("CQTheta");

    // const tensorField& CQDTheta =
    //     patch().lookupPatchField<surfaceTensorField, tensor>("CQDTheta");

    const vectorField explicitQ =
        patch().lookupPatchField<surfaceVectorField, vector>("explicitQ");

    const scalarField delta(1.0/patch().deltaCoeffs());

    // This is the old coefficient in force boundary condition
    // whose inverse was being multiplied
    const tensorField A(CQW/delta);
    
    // This coefficent changes due to the spring
    const tensorField invAmodified
    (
        inv
        (
            CQW/delta - Ks_
        )
    );

    // Modified displacement at boundary due to spring
    const vectorField newDW
    (
     (invAmodified & (force() - explicitQ))
      - (invAmodified & (CQTheta & DTheta))
      + (invAmodified & (A & DW.patchInternalField()))
    );
    
    DW = newDW;

    // Setting W_ field
    fixedValueFvPatchField<vector>::operator==((*this) + newDW);

    fixedValueFvPatchVectorField::evaluate();
}


// Write
void linearSpringForceBeamDisplacementNRFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);

    os.writeKeyword("linearSpringCoeffs") << nl;
    os << token::BEGIN_BLOCK << nl;
    os.writeKeyword("springStiffness") << springStiffness_ << token::END_STATEMENT << nl;
    os.writeKeyword("springDirection") << springDirection_ << token::END_STATEMENT << nl;
    os.writeKeyword("springForce") << force() << token::END_STATEMENT << nl;
    os << token::END_BLOCK << nl;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    linearSpringForceBeamDisplacementNRFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
