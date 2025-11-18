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
    gDamping_
    (
        readScalar(beamMomentumContribDict_.lookup("gDamping"))
    ),
    gStiffness_
    (
        readScalar(beamMomentumContribDict_.lookup("gStiffness"))
    ),
    groundZ_
    (
        readScalar(beamMomentumContribDict_.lookup("groundZ"))
    )
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
    // Field<scalarSquareMatrix>& result = tresult.ref();

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
                dimensionedVector("refW", dimVelocity, vector::zero)
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
    // ALSO need to add friction part of contact
    label cellsInContact = 0;

    forAll(result, cellI)
    {
        const vector coord = refW[cellI] + W[cellI];
        Info<< "coord " << coord << endl;
        if (coord.z() < groundZ_)
        {
            cellsInContact++;
            result[cellI][vector::Z] +=
                (2.0*gStiffness_*R*(coord.z() - groundZ_))
              - (2.0*gDamping_*R*max(U[cellI].component(2), 0));
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
