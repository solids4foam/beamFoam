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

#include "groundContactContribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(groundContactContribution, 0);
addToRunTimeSelectionTable
(
    beamMomentumContribution, groundContactContribution, dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


groundContactContribution::groundContactContribution
(
    const word& name,
    const dictionary& dict
)
:
    beamMomentumContribution
    (
        name,
        dict
    ),
    beamMomentumContribDict_(dict.subDict(name + "Coeffs")),
    kNormal_
    (
        readScalar(beamMomentumContribDict_.lookup("kNormal"))
    ),
    cNormal_
    (
        readScalar(beamMomentumContribDict_.lookup("cNormal"))
    ),
    kTangential_
    (
        readScalar(beamMomentumContribDict_.lookup("kTangential"))
    ),
    muFriction_
    (
        readScalar(beamMomentumContribDict_.lookup("muFriction"))
    ),
    groundZ_
    (
        readScalar(beamMomentumContribDict_.lookup("groundZ"))
    ),
    contactStiffness_()
{
    Info<< "Found beamMomentumContribution type: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<Field<scalarSquareMatrix>> groundContactContribution::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<Field<scalarSquareMatrix>> tresult
    (
        new Field<scalarSquareMatrix>
        (
            mesh.nCells(), scalarSquareMatrix(6, 0.0)
        )
    );
    Field<scalarSquareMatrix>& result = tresult.ref();

    // The diagonal holds the derivative of the residual (internal plus
    // external force minus inertia) with respect to W, like the implicit
    // inertia coefficient. The normal contact force is
    // 2*R*kNormal*(groundZ - z)*L, so its derivative with respect to W.z is
    // -2*R*kNormal*L. The damping and friction terms stay explicit.
    // contactStiffness_ is set by linearMomentumSource, which the beam model
    // calls first in the same Newton iteration
    if (contactStiffness_.size() == mesh.nCells())
    {
        forAll(result, cellI)
        {
            result[cellI](2, 2) -= contactStiffness_[cellI];
        }
    }

    return tresult;
}


tmp<vectorField> groundContactContribution::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    IOobject refWHeader
    (
        "refW",
        "0",
        mesh,
        IOobject::MUST_READ
    );

    autoPtr<volVectorField> refWPtr;

    if (refWHeader.typeHeaderOk<volVectorField>(true))
    {
        // Info<< "Reading refW from 0/" << endl;

        refWPtr.reset
        (
            new volVectorField(refWHeader, mesh)
        );
    }
    else
    {
        // Info<< "refW not found → using default zero field" << endl;

        refWPtr.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "refW",
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedVector("refW", dimLength, vector::zero)
            )
        );
    }

    volVectorField& refW = refWPtr();

    // Beam Radius
    const scalar R = bm.R();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    vectorField& result = tresult.ref();

    // TO-DO: Make the ground contact not just for z-direction but any
    // user specified direction

    // Create spline using current beam points and tangents data
    HermiteSpline spline
    (
        bm.currentBeamPoints(),
        bm.currentBeamTangents()
    );

    // Evaluate dRdS - tangents to beam centreline at beam CV cell-centres
    const vectorField& dRdScell = spline.midPointDerivatives();

    // Tangential component of velocity vector
    vectorField Ut
    (
        (
            (U.internalField() & dRdScell)
            *dRdScell
        )
    );

    vectorField UtHat (Ut/(mag(Ut) + SMALL));

    // Beam cell lengths: the forces below are per unit length
    const volScalarField& L = bm.L();

    contactStiffness_.setSize(mesh.nCells());
    contactStiffness_ = 0;

    // Initialise beam cells in contact with ground
    label cellsInContact = 0;

    forAll(result, cellI)
    {
        const vector coord = refW[cellI] + W[cellI];
        // Info<< "coord " << coord << endl;
        if (coord.z() < groundZ_)
        {
            cellsInContact++;
            contactStiffness_[cellI] = 2*R*kNormal_*L[cellI];

            // Spring plus damper, damping motion in both directions; the
            // ground can push but not pull
            const scalar f_gc_normal = max
            (
                (2*R*kNormal_*(groundZ_ - coord.z()))
              - (2*R*cNormal_*U[cellI].component(2)),
                0
            );

            vector f_gc_tangential(vector::zero);

            if
            (
                kTangential_*2.0*R*mag(Ut[cellI])
             >= muFriction_*mag(f_gc_normal)
            )
            {
                // Max value of friction as per Coulomb's law
                f_gc_tangential =
                    -muFriction_*mag(f_gc_normal)*UtHat[cellI];
            }
            else
            {
                f_gc_tangential = -2.0*R*kTangential_*Ut[cellI];
            }
            // The beam model adds this to its source, which holds minus the
            // external force on each cell (as for gravity and the
            // distributed load q), so store minus the force times the
            // cell length
            result[cellI][vector::X] -= f_gc_tangential.x()*L[cellI];
            result[cellI][vector::Y] -= f_gc_tangential.y()*L[cellI];
            result[cellI][vector::Z] -=
                (f_gc_normal + f_gc_tangential.z())*L[cellI];

        }
     }

    Info<< "Number of cells in contact : " << cellsInContact << endl;

    return tresult;
}


tmp<vectorField> groundContactContribution::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    // vectorField& result = tresult.ref();

    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
