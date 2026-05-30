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

#include "almDragContribution.H"
#include "addToRunTimeSelectionTable.H"
#include "HermiteSpline.H"
#include "samplingFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(almDragContribution, 0);
addToRunTimeSelectionTable
(
    beamMomentumContribution, almDragContribution, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

almDragContribution::almDragContribution
(
    const word& name,
    const dictionary& dict
)
:
    beamMomentumContribution(name, dict),
    coeffs_(dict.subDict(name + "Coeffs")),
    Cdn_(readScalar(coeffs_.lookup("Cdn"))),
    Cdt_(readScalar(coeffs_.lookup("Cdt"))),
    forceRelaxation_
    (
        coeffs_.lookupOrDefault<scalar>("forceRelaxation", 1.0)
    ),
    samplingRadius_
    (
        coeffs_.lookupOrDefault<scalar>("samplingRadius", 5.0)
    ),
    samplingReferenceLength_
    (
        coeffs_.lookupOrDefault<scalar>("almSamplingReferenceLength", -1.0)
    ),
    dragReferenceLength_
    (
        coeffs_.lookupOrDefault<scalar>("almDragReferenceLength", -1.0)
    ),
    samplingPlaneNormal_(coeffs_.lookup("almSamplingPlaneNormal")),
    upstreamDir_(coeffs_.lookup("almUpstreamDir")),
    groundContactActive_
    (
        coeffs_.lookupOrDefault<bool>("groundContactActive", false)
    ),
    groundZ_
    (
        coeffs_.lookupOrDefault<scalar>("groundZ", 0.0)
    ),
    searchEnginePtr_()
{
    Info<< "Found beamMomentumContribution type: " << typeName << endl;

    if (mag(samplingPlaneNormal_) <= SMALL)
    {
        FatalIOErrorInFunction(coeffs_)
            << "almSamplingPlaneNormal magnitude must be positive. "
            << "Current value: " << samplingPlaneNormal_
            << exit(FatalIOError);
    }

    if (mag(upstreamDir_) <= SMALL)
    {
        FatalIOErrorInFunction(coeffs_)
            << "almUpstreamDir magnitude must be positive. "
            << "Current value: " << upstreamDir_
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void almDragContribution::preEvolve(const beamModel& bm)
{
    const fvMesh& beamMesh = bm.solutionW().mesh();

    const fvMesh& fluidMesh =
        beamMesh.time().db().parent().lookupObject<fvMesh>("region0");

    if (!searchEnginePtr_.valid())
    {
        searchEnginePtr_.reset(new meshSearch(fluidMesh));
    }

    // The "almForce" field is owned by the beam model and registered on the
    // beam mesh.  Look it up by name to avoid coupling to a concrete derived
    // beam type.
    volVectorField& almForce =
        beamMesh.lookupObjectRef<volVectorField>("almForce");

    // Per-cell line force at the start of this evolve() - used for
    // outer-corrector under-relaxation
    const vectorField almForceStart(almForce.internalField());

    // Beam-aligned tangent at beam cell centres (broadcast from master)
    vectorField dRdScell;

    if (Pstream::master())
    {
        HermiteSpline spline
        (
            bm.currentBeamPoints(),
            bm.currentBeamTangents()
        );

        dRdScell = spline.midPointDerivatives();
    }

    label nDRdS = dRdScell.size();
    Pstream::broadcast(nDRdS);

    if (!Pstream::master())
    {
        dRdScell.setSize(nDRdS, vector::zero);
    }

    Pstream::broadcast(dRdScell);

    const vectorField dRdScellHat(dRdScell/(mag(dRdScell) + SMALL));

    // Beam cell-centre global coordinates (reference + total displacement)
    const volVectorField& W = bm.solutionW();
    const vectorField beamCellCenterCoord =
        beamMesh.C() + W.mesh().lookupObject<volVectorField>("refW") + W;

    labelList seedCellIDs(beamMesh.nCells(), -1);

    const scalar samplingRefL =
        (samplingReferenceLength_ > 0) ? samplingReferenceLength_ : bm.R();

    if (samplingRefL <= SMALL)
    {
        FatalErrorInFunction
            << "almSamplingReferenceLength must be positive. "
            << "Current value: " << samplingRefL
            << abort(FatalError);
    }

    const scalar dragRefD =
        (dragReferenceLength_ > 0) ? dragReferenceLength_ : 2.0*bm.R();

    if (dragRefD <= SMALL)
    {
        FatalErrorInFunction
            << "almDragReferenceLength must be positive. "
            << "Current value: " << dragRefD
            << abort(FatalError);
    }

    std::tuple
    <
        tmp<volVectorField>,
        tmp<volScalarField>,
        tmp<volVectorField>,
        labelList
    > fluidInfo =
        getFluidVelocity
        (
            fluidMesh,
            beamMesh,
            beamCellCenterCoord,
            seedCellIDs,
            groundZ_,
            groundContactActive_,
            dRdScell,
            samplingRefL,
            samplingRadius_,
            samplingPlaneNormal_,
            upstreamDir_,
            searchEnginePtr_()
        );

    const volVectorField fluidU(std::get<0>(fluidInfo).ref());
    const volScalarField cellMarker(std::get<1>(fluidInfo).ref());

    if (beamMesh.time().writeTime())
    {
        cellMarker.write();
    }

    // Relative velocity (fluid minus beam)
    const volVectorField& U =
        beamMesh.lookupObject<volVectorField>("U");
    const vectorField Urel(fluidU.internalField() - U.internalField());

    const vectorField Ut
    (
        (Urel & dRdScellHat)*dRdScellHat
    );

    const vectorField Un
    (
        Urel - Ut
    );

    const scalar rhoF = bm.rhoFluid().value();
    const scalar coeffN = 0.5*rhoF*Cdn_*dragRefD;
    const scalar coeffT = 0.5*rhoF*Cdt_*dragRefD;

    vectorField& almForceI = almForce.primitiveFieldRef();

    forAll(almForceI, cellI)
    {
        const vector un = Un[cellI];
        const vector ut = Ut[cellI];
        const vector vs = un + ut;
        const scalar vsMag = mag(vs);

        const vector fnLineTarget = -coeffN*vsMag*un;
        const vector ftLineTarget = -coeffT*vsMag*ut;
        const vector FLineTarget = fnLineTarget + ftLineTarget;

        almForceI[cellI] =
            forceRelaxation_*FLineTarget
          + (1.0 - forceRelaxation_)*almForceStart[cellI];
    }

    almForce.correctBoundaryConditions();

    const vector Fsum(sum(almForce*bm.L()).value());
    Info<< "sum F on the beam (relaxed) = "
        << "(" << Fsum.x() << ", "
        << Fsum.y() << ", "
        << Fsum.z() << ") N" << nl;
}


tmp<Field<scalarSquareMatrix>> almDragContribution::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const fvMesh& mesh = bm.solutionW().mesh();

    return tmp<Field<scalarSquareMatrix>>
    (
        new Field<scalarSquareMatrix>
        (
            mesh.nCells(), scalarSquareMatrix(6, 0.0)
        )
    );
}


tmp<vectorField> almDragContribution::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const fvMesh& beamMesh = bm.solutionW().mesh();

    const volVectorField& almForce =
        beamMesh.lookupObject<volVectorField>("almForce");

    const volScalarField& L = bm.L();

    tmp<vectorField> tresult(new vectorField(beamMesh.nCells(), vector::zero));
    vectorField& result = tresult.ref();

    const vectorField& almForceI = almForce.internalField();
    const scalarField& LI = L.internalField();

    forAll(result, cellI)
    {
        result[cellI] = almForceI[cellI]*LI[cellI];
    }

    return tresult;
}


tmp<vectorField> almDragContribution::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const fvMesh& mesh = bm.solutionW().mesh();

    return tmp<vectorField>
    (
        new vectorField(mesh.nCells(), vector::zero)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
