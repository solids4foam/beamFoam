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

#include "morisonForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(morisonForce, 0);
defineRunTimeSelectionTable(morisonForce, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<morisonForce> morisonForce::New
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& solidDict
)
{
    word morisonForceTypeName;

    // Enclose the creation of the waveProp to ensure it is
    // deleted before the morisonForce is created otherwise the dictionary
    // is entered in the database twice
    {
        if (runTime.db().foundObject<IOdictionary>("forceProperties"))
        {
            runTime.db().lookupObject<IOdictionary>("forceProperties").lookup("forceType")
                 >> morisonForceTypeName;
        }
        else
        {
            IOdictionary dict
            (
                IOobject
                (
                    "forceProperties",
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            dict.lookup("forceType") >> morisonForceTypeName;
        }
    }

    auto* ctorPtr = dictionaryConstructorTable(morisonForceTypeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            solidDict,
            "morisonForce",
            morisonForceTypeName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<morisonForce>(ctorPtr(runTime, mesh, solidDict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


morisonForce::morisonForce
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& solidDict
)
:
    IOdictionary
    (
        IOobject
        (
            "forceProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    runTime_(runTime),

    mesh_(mesh),

    dragCoeff_
    (
        this->subDict("forceCoeffs").lookup("CD")
    ),
    inertiaCoeff_
    (
        this->subDict("forceCoeffs").lookup("CM")
    ),
    frictionCoeff_
    (
        this->subDict("forceCoeffs").lookup("Cf")
    ),
    rhoFluid_
    (
        this->subDict("forceCoeffs").lookup("rhoFluid")
    ),
    rhoStem_
    (
        solidDict.lookup("rho")
    ),
    G_(9.81),
    vertDir_
    (
        this->subDict("forceCoeffs").lookup("posVertDir")
    ),

    addedMassCoeffs_
    (
        IOobject
        (
            "addedMass",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("null", dimMass, symmTensor::zero),
        "zeroGradient"
    ),

    dragCoeffs_
    (
        IOobject
        (
            "drag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("null", dimMass/dimTime, symmTensor::zero),
        "zeroGradient"
    ),

    linearize_(this->subDict("forceCoeffs").lookupOrDefault<Switch>("linearize", false))
{
}


morisonForce::~morisonForce()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void morisonForce::updateForce
(
    volVectorField& q,
    const volVectorField& dW,
    const volVectorField& U,
    const tmp<vectorField> tangentFace,
    const scalar& width,
    const scalar& area,
    const volScalarField& L
) const
{
    // Initial stuff
    const vectorField& cc(mesh_.cellCentres());
    scalar t = runTime_.time().value();
    scalar dt =  runTime_.deltaT().value();

    // Make the explicit acceleration calculation
    volVectorField stemAcc(fvc::ddt(U));

    // Reference to force field
    vectorField& qi(q.primitiveFieldRef());

    // Make the tanget field
    vectorField tangentCell(qi.size(), vector::zero);
    vectorField normalCell(qi.size(), vector::zero);

    forAll (tangentCell, celli)
    {
        tangentCell[celli] = 0.5*(tangentFace()[celli] + tangentFace()[celli + 1]);
        tangentCell[celli] /= Foam::mag(tangentCell[celli]);

        // Hard-coded (again) with forcing in y-direction
        if (linearize_)
        {
            normalCell[celli] = vector(1,0,0);
        }
        else
        {
            normalCell[celli] = vector(0, 1, 0) ^ tangentCell[celli];
            normalCell[celli] /= Foam::mag(normalCell[celli]);
        }
    }

    // Initial stuff
    scalar factDrag(0.5*rhoFluid_.value()*dragCoeff_.value()*width);
    scalar factIner
    (
        rhoFluid_.value()*inertiaCoeff_.value()*Foam::sqr(width)/4.0*M_PI
    );
    scalar factFroudeKrylov(rhoFluid_.value()*area);

    // There are two sides
    scalar factFric(2*frictionCoeff_.value()*rhoFluid_.value()*width);

    // Evaluate the added mass and inertia
    symmTensorField& am(addedMassCoeffs_.primitiveFieldRef());
    symmTensorField& dc(dragCoeffs_.primitiveFieldRef());

    const scalarField& l(L.primitiveField());

    vectorField explicitIner(am.size(), vector::zero);
    vectorField explicitDrag(dc.size(), vector::zero);

    forAll (am, celli)
    {
        am[celli] = factIner*l[celli]*Foam::sqr(normalCell[celli]);
        dc[celli] = factDrag*l[celli]*Foam::sqr(normalCell[celli]);

        // Add all contributions to the explicit contributions
        vector& ei(explicitIner[celli]);
        vector& ed(explicitDrag[celli]);

        ei.x() = Foam::neg(am[celli].xy())*am[celli].xy()
            + Foam::neg(am[celli].xz())*am[celli].xz();
        ei.y() = Foam::neg(am[celli].xy())*am[celli].xy()
            + Foam::neg(am[celli].yz())*am[celli].yz();
        ei.z() = Foam::neg(am[celli].xz())*am[celli].xz()
            + Foam::neg(am[celli].yz())*am[celli].yz();

        ed.x() = Foam::neg(dc[celli].xy())*dc[celli].xy()
            + Foam::neg(dc[celli].xz())*dc[celli].xz();
        ed.y() = Foam::neg(dc[celli].xy())*dc[celli].xy()
            + Foam::neg(dc[celli].yz())*dc[celli].yz();
        ed.z() = Foam::neg(dc[celli].xz())*dc[celli].xz()
            + Foam::neg(dc[celli].yz())*dc[celli].yz();

        // Remove all the negative contributions from inertia
        am[celli].xy() *= Foam::pos(am[celli].xy());
        am[celli].xz() *= Foam::pos(am[celli].xz());
        am[celli].yz() *= Foam::pos(am[celli].yz());

        // Remove all the negative contributions from drag
        dc[celli].xy() *= Foam::pos(dc[celli].xy());
        dc[celli].xz() *= Foam::pos(dc[celli].xz());
        dc[celli].yz() *= Foam::pos(dc[celli].yz());
    }

    // Loop over all cells for the explicit terms and the ||u - phi,t||_2
    // multiplier on the explicit drag force
//    point XX(point(0.0, 0.0, 0.3));
//
//    vector tmpU = this->U(XX, t).value();
//    vector tmpA = this->acc(XX, t).value();
//
//    Info << "U/A:" << tab << t << tab << tmpU.x() << tab << tmpU.z() << tab <<
//            tmpA.x() << tab << tmpA.z() << endl;

    forAll (cc, celli)
    {
        dimensionedVector u (this->U(cc[celli], t));
        dimensionedVector acc(this->acc(cc[celli], t));

        scalar relu = Foam::mag((u.value() - U.internalField()[celli]) & normalCell[celli]);
//        scalar relu = Foam::mag((u.value() - U.internalField()[celli]));

        vector drag = vector::zero;
        vector utau((u.value() & tangentCell[celli])*tangentCell[celli]);

        u.value()   = (u.value() & normalCell[celli])*normalCell[celli];
        acc.value() = (acc.value() & normalCell[celli])*normalCell[celli];

        if (linearize_)
        {
            drag = factDrag*this->Uabsolute()*u.value();
            dc[celli] *= this->Uabsolute();
            explicitDrag[celli] *= this->Uabsolute();
        }
        else
        {
            drag = factDrag*relu*u.value();
            dc[celli] *= relu;
            explicitDrag[celli] *= relu;
        }

        vector inertia = (factIner + factFroudeKrylov)*acc.value();
        vector friction = factFric*Foam::mag(utau)*utau;

        // Add the explicit drag term to the formulation
        drag -= Foam::cmptMultiply(explicitDrag[celli], U.internalField()[celli]);
        inertia -= Foam::cmptMultiply(explicitIner[celli], stemAcc.internalField()[celli]);

        qi[celli] = drag + inertia + friction
            + (area*G_*(rhoFluid_.value() - rhoStem_.value()))*vertDir_;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
