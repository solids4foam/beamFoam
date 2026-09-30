/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

Description
    Coupling between the beam and a rigidBodyEnd owned by beamFoam.

    rigidBodyCoupling selects the method:
        none         no rigid body (default without a rigidBodyEnd dict)
        partitioned  body and beam solved in turn, with optional coupling
                     iterations (default with a rigidBodyEnd dict)
        monolithic   body solved inside the BlockEigen system (not yet
                     implemented: Phase 1 of codexLogs/monolithicCouplingPlan)

    The legacy blockEigen* switches are used by the moorFV-driven path and
    cannot be combined with rigidBodyEnd.

\*---------------------------------------------------------------------------*/

#include "coupledTotalLagNewtonRaphsonBeam.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledTotalLagNewtonRaphsonBeam::readRigidBodyEnd()
{
    const dictionary& coeffs = beamProperties();
    const bool hasRigidBodyEnd = coeffs.found("rigidBodyEnd");

    rigidBodyCoupling_ =
        coeffs.getOrDefault<word>
        (
            "rigidBodyCoupling",
            hasRigidBodyEnd ? "partitioned" : "none"
        );

    if
    (
        !List<word>({"none", "partitioned", "monolithic"})
            .found(rigidBodyCoupling_)
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "Invalid rigidBodyCoupling " << rigidBodyCoupling_ << nl
            << "Valid options are none, partitioned and monolithic"
            << exit(FatalIOError);
    }

    const wordList legacySwitches
    ({
        "blockEigenForceCoupling",
        "blockEigenKinematicCoupling"
    });

    bool legacyCouplingOn = false;
    forAll(legacySwitches, i)
    {
        legacyCouplingOn =
            legacyCouplingOn
         || coeffs.getOrDefault<Switch>(legacySwitches[i], false);
    }

    if (!hasRigidBodyEnd)
    {
        if (rigidBodyCoupling_ != "none")
        {
            FatalIOErrorInFunction(coeffs)
                << "rigidBodyCoupling " << rigidBodyCoupling_
                << " needs a rigidBodyEnd dictionary"
                << exit(FatalIOError);
        }

        if (legacyCouplingOn)
        {
            Info<< "Legacy blockEigen coupling switches are set. They are "
                << "used by the moorFV-driven path and do not give a "
                << "monolithic beam / rigid-body solve." << endl;
        }

        return;
    }

    if (rigidBodyCoupling_ == "none")
    {
        FatalIOErrorInFunction(coeffs)
            << "A rigidBodyEnd dictionary is given but rigidBodyCoupling "
            << "is none" << exit(FatalIOError);
    }

    if (rigidBodyCoupling_ == "monolithic")
    {
        FatalIOErrorInFunction(coeffs)
            << "rigidBodyCoupling monolithic is not implemented yet "
            << "(Phase 1 of codexLogs/monolithicCouplingPlan.pdf). "
            << "Use partitioned." << exit(FatalIOError);
    }

    if (legacyCouplingOn)
    {
        FatalIOErrorInFunction(coeffs)
            << "rigidBodyEnd cannot be combined with the legacy "
            << legacySwitches << " switches" << exit(FatalIOError);
    }

    const word defaultPatchName =
        endPatchIndex() >= 0
      ? mesh().boundaryMesh()[endPatchIndex()].name()
      : word("right");

    rigidBodyEndPtr_.reset
    (
        new rigidBodyEnd
        (
            runTime(),
            coeffs.subDict("rigidBodyEnd"),
            defaultPatchName,
            g().value()
        )
    );

    const label patchI =
        mesh().boundaryMesh().findPatchID(rigidBodyEndPtr_().patchName());

    if (patchI < 0)
    {
        FatalErrorInFunction
            << "rigidBodyEnd patch " << rigidBodyEndPtr_().patchName()
            << " not found" << abort(FatalError);
    }

    if (!isA<fixedValueFvPatchVectorField>(W_.boundaryField()[patchI]))
    {
        FatalErrorInFunction
            << "rigidBodyEnd needs a fixedValue W boundary condition on "
            << "patch " << rigidBodyEndPtr_().patchName()
            << abort(FatalError);
    }

    if (W_.boundaryField()[patchI].size() != 1)
    {
        FatalErrorInFunction
            << "rigidBodyEnd patch " << rigidBodyEndPtr_().patchName()
            << " must have exactly one face on this processor"
            << abort(FatalError);
    }

    rigidBodyEndPtr_().setAttachmentOffset(W_.boundaryField()[patchI][0]);

    Info<< "rigidBodyCoupling " << rigidBodyCoupling_ << endl;
}


scalar coupledTotalLagNewtonRaphsonBeam::evolve()
{
    if (rigidBodyEndPtr_.valid())
    {
        if (rigidBodyDataValid_)
        {
            FatalErrorInFunction
                << "rigidBodyEnd cannot be used when rigid-body data is also "
                << "supplied by moorFV" << abort(FatalError);
        }

        return evolvePartitionedRigidBodyEnd();
    }

    return evolveBeam();
}


scalar coupledTotalLagNewtonRaphsonBeam::evolvePartitionedRigidBodyEnd()
{
    rigidBodyEnd& body = rigidBodyEndPtr_();

    body.newTimeStep();

    const label patchI = mesh().boundaryMesh().findPatchID(body.patchName());

    // Coupling starts from the beam force of the last accepted state
    vector beamForce = body.beamForce();
    vector x = body.displacement();

    scalar initialResidual = 0;
    label nBeamSolves = 0;

    for (label iter = 1; iter <= body.nCouplingIterations(); ++iter)
    {
        const vector xNewmark = body.newmarkDisplacement(beamForce);

        // Change relative to the larger of the step and total displacement,
        // so the test still works when the body is at rest
        const scalar change =
            mag(xNewmark - x)
           /max
            (
                max(mag(xNewmark - body.oldDisplacement()), mag(xNewmark)),
                VSMALL
            );

        if (iter > 1)
        {
            Info<< "rigidBodyEnd coupling iteration " << iter - 1
                << ": relative displacement change " << change << endl;

            if (change < body.couplingTolerance())
            {
                break;
            }
        }

        x =
            iter > 1
          ? body.couplingRelaxation()*xNewmark
          + (1 - body.couplingRelaxation())*x
          : xNewmark;

        W_.boundaryFieldRef()[patchI] == (x + body.attachmentOffset());

        const scalar residual = evolveBeam();

        if (nBeamSolves == 0)
        {
            initialResidual = residual;
        }

        ++nBeamSolves;

        // Force from the beam on the body
        beamForce = -Q_.boundaryField()[patchI][0];
    }

    body.accept(x, beamForce, nBeamSolves);
    body.write();

    Info<< "rigidBodyEnd: displacement " << body.displacement()
        << ", velocity " << body.velocity()
        << ", beam force " << body.beamForce()
        << ", beam solves " << nBeamSolves << endl;

    return initialResidual;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
