/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Partitioned coupling between the beam and a rigidBodyEnd: the body and
    the beam are solved in turn, with coupling iterations within each time
    step (see rigidBodyEnd.H).

    The body sets the beam end displacement through a fixedValue W
    condition on its patch (and, for a welded joint, the end rotation
    through a fixedValue Theta condition). The beam returns the end-face
    force (and, for a welded joint, the end moment).

\*---------------------------------------------------------------------------*/

#include "coupledTotalLagNewtonRaphsonBeam.H"
#include "fixedValueFvPatchFields.H"
#include "momentBeamRotationNRFvPatchVectorField.H"
#include "pseudoVector.H"

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

    if (!coeffs.found("rigidBodyEnd"))
    {
        return;
    }

    if (Pstream::parRun())
    {
        FatalIOErrorInFunction(coeffs)
            << "rigidBodyEnd does not support parallel runs yet: run the "
            << "case in serial" << exit(FatalIOError);
    }

    const word coupling
    (
        coeffs.getOrDefault<word>("rigidBodyCoupling", "partitioned")
    );

    if (coupling != "partitioned")
    {
        FatalIOErrorInFunction(coeffs)
            << "Invalid rigidBodyCoupling " << coupling << nl
            << "Valid options are: partitioned" << exit(FatalIOError);
    }

    // The body is advanced with the beam's Newmark parameters, which are
    // only set for the Newmark scheme
    if (d2dt2SchemeName_ != "Newmark")
    {
        FatalIOErrorInFunction(coeffs)
            << "rigidBodyEnd needs the Newmark d2dt2Scheme in system/fvSchemes"
            << ", got " << d2dt2SchemeName_ << exit(FatalIOError);
    }

    const word defaultPatchName =
        endPatchIndex() >= 0
      ? mesh().boundaryMesh()[endPatchIndex()].name()
      : word("right");

    rigidBodyEndPtr_.reset
    (
        new rigidBodyEnd
        (
            mesh(),
            coeffs.subDict("rigidBodyEnd"),
            defaultPatchName,
            g().value(),
            betaN_,
            gammaN_
        )
    );

    const rigidBodyEnd& body = rigidBodyEndPtr_();

    const label patchI = mesh().boundaryMesh().findPatchID(body.patchName());

    if (patchI < 0)
    {
        FatalErrorInFunction
            << "rigidBodyEnd patch " << body.patchName() << " not found"
            << abort(FatalError);
    }

    if (W_.boundaryField()[patchI].size() != 1)
    {
        FatalErrorInFunction
            << "rigidBodyEnd patch " << body.patchName() << " has "
            << W_.boundaryField()[patchI].size() << " faces: it must be "
            << "the end of a single beam (one face)" << abort(FatalError);
    }

    // The body must be the only thing that sets W on its patch. Conditions
    // derived from fixedValue (e.g. fixedBeamDisplacement, the force
    // conditions) reset or recompute the value, so the exact type is
    // required
    const word& WType = W_.boundaryField()[patchI].type();

    if (WType != fixedValueFvPatchVectorField::typeName)
    {
        FatalErrorInFunction
            << "rigidBodyEnd needs a fixedValue W boundary condition on "
            << "patch " << body.patchName() << ", got " << WType
            << abort(FatalError);
    }

    const word& ThetaType = Theta_.boundaryField()[patchI].type();

    if (body.welded())
    {
        if (ThetaType != fixedValueFvPatchVectorField::typeName)
        {
            FatalErrorInFunction
                << "A welded rigidBodyEnd needs a fixedValue Theta boundary "
                << "condition on patch " << body.patchName() << ", got "
                << ThetaType << abort(FatalError);
        }
    }
    else
    {
        // A pin passes no moment to the body, so the beam end must be free
        // to rotate: a zero moment condition. A condition that fixes the
        // rotation would clamp the beam end, and the end moment it carries
        // would never reach the body
        if (ThetaType != momentBeamRotationNRFvPatchVectorField::typeName)
        {
            FatalErrorInFunction
                << "A pinned rigidBodyEnd needs a momentBeamRotationNR Theta "
                << "boundary condition with zero moment on patch "
                << body.patchName() << ", so that the beam end is free to "
                << "rotate; got " << ThetaType << abort(FatalError);
        }

    }

    rigidBodyEndPatchIndex_ = patchI;

    if (!body.welded())
    {
        checkRigidBodyEndPinMoment();
    }

    // The end position without W includes the reference displacement refWf
    // (e.g. the placement from setInitialPositionBeam). The beam force and
    // moment on the body are minus those on the end face
    rigidBodyEndPtr_().attach
    (
        mesh().Cf().boundaryField()[patchI][0]
      + refWf_.boundaryField()[patchI][0],
        W_.boundaryField()[patchI][0],
        -Q_.boundaryField()[patchI][0],
        rigidBodyEndBeamMoment()
    );

    Info<< "rigidBodyCoupling " << coupling << endl;
}


scalar coupledTotalLagNewtonRaphsonBeam::evolvePartitionedRigidBodyEnd()
{
    rigidBodyEnd& body = rigidBodyEndPtr_();

    body.newTimeStep();
    if (body.welded())
    {
        storeWeldLambda0();
    }

    const label patchI = rigidBodyEndPatchIndex_;

    // Length that converts a rotation change to a displacement change. A
    // welded body can turn without its arm moving (e.g. a centre of mass
    // at the attachment), so its radius of gyration sets a lower bound
    const scalar armLength =
        body.welded()
      ? max(mag(body.arm(vector::zero)), body.radiusOfGyration())
      : mag(body.arm(vector::zero));

    // Coupling starts from the beam force and moment of the last accepted
    // state
    vector beamForce = body.beamForce();
    vector beamMoment = body.beamMoment();
    vector x = body.displacement();
    vector theta = vector::zero;

    scalar initialResidual = 0;
    label nBeamSolves = 0;
    bool converged = false;
    scalar change = GREAT;

    // Each pass moves the body with the latest beam force and moment, then
    // solves the beam with the end set by the body. The pass after the last
    // beam solve only measures the change, to tell whether it converged
    for (label iter = 1; ; ++iter)
    {
        const vector xNewmark = body.newmarkDisplacement(beamForce);
        const vector thetaNewmark =
            body.newmarkRotation(beamForce, theta, beamMoment);

        if (iter > 1)
        {
            // Change relative to the larger of the step and total motion,
            // so the test still works when the body is at rest. A welded
            // body's total motion includes its rotation, which may be all
            // it does
            const scalar totalMotion =
                mag(xNewmark)
              + (
                    body.welded()
                  ? armLength*body.totalRotation(thetaNewmark)
                  : 0
                );

            change =
                (mag(xNewmark - x) + armLength*mag(thetaNewmark - theta))
               /max
                (
                    max
                    (
                        mag(xNewmark - body.oldDisplacement())
                      + armLength*mag(thetaNewmark),
                        totalMotion
                    ),
                    VSMALL
                );

            Info<< "rigidBodyEnd coupling iteration " << iter - 1
                << ": relative change " << change << endl;

            if (change < body.couplingTolerance())
            {
                converged = true;
                break;
            }

            if (nBeamSolves >= body.nCouplingIterations())
            {
                break;
            }
        }

        const scalar w = iter > 1 ? body.couplingRelaxation() : 1;
        x = w*xNewmark + (1 - w)*x;
        theta = w*thetaNewmark + (1 - w)*theta;

        W_.boundaryFieldRef()[patchI][0] = rigidBodyEndAttachmentW(x, theta);

        if (body.welded())
        {
            setWeldRotationIncrement
            (
                weldRotationIncrement(body.stepRotation(theta))
            );
        }

        const scalar residual = evolveBeam();

        if (nBeamSolves == 0)
        {
            initialResidual = residual;
        }

        ++nBeamSolves;

        if (!body.welded())
        {
            checkRigidBodyEndPinMoment();
        }

        // Force and moment from the beam on the body
        beamForce = -Q_.boundaryField()[patchI][0];
        beamMoment = rigidBodyEndBeamMoment();
    }

    if (!converged)
    {
        WarningInFunction
            << "rigidBodyEnd coupling did not converge in "
            << nBeamSolves << " beam solves at time "
            << runTime().timeName() << ": relative change " << change
            << ", couplingTolerance " << body.couplingTolerance() << nl
            << "    Increase nCouplingIterations, or reduce deltaT or use "
            << "couplingRelaxation if the change does not decrease (see "
            << "rigidBodyEnd.H)" << endl;
    }

    body.accept(x, theta, beamForce, beamMoment, nBeamSolves, converged);
    body.write();

    Info<< "rigidBodyEnd: displacement " << body.displacement()
        << ", velocity " << body.velocity()
        << ", angular velocity " << body.angularVelocity()
        << ", beam force " << body.beamForce()
        << ", beam moment " << body.beamMoment()
        << ", beam solves " << nBeamSolves
        << (converged ? "" : " (not converged)") << endl;

    return initialResidual;
}


vector coupledTotalLagNewtonRaphsonBeam::rigidBodyEndAttachmentW
(
    const vector& x,
    const vector& theta
) const
{
    const rigidBodyEnd& body = rigidBodyEndPtr_();

    return body.attachmentDisplacement(x, theta) + body.attachmentOffset();
}


void coupledTotalLagNewtonRaphsonBeam::storeWeldLambda0()
{
    if (runTime().timeIndex() == weldTimeIndex_)
    {
        return;
    }

    weldTimeIndex_ = runTime().timeIndex();
    weldLambda0_ = Lambdaf_.boundaryField()[rigidBodyEndPatchIndex_][0];
}


vector coupledTotalLagNewtonRaphsonBeam::weldRotationIncrement
(
    const tensor& R
) const
{
    // Lambdaf_ excludes the reference rotation refLambdaf_, which is the
    // same on both sides and cancels. The increment is one coupling
    // iteration's turn of the end face, so always small
    return pseudoVector
    (
        R
      & weldLambda0_
      & Lambdaf_.boundaryField()[rigidBodyEndPatchIndex_][0].T()
    );
}


void coupledTotalLagNewtonRaphsonBeam::setWeldRotationIncrement
(
    const vector& increment
)
{
    const label patchI = rigidBodyEndPatchIndex_;

    // The fixedValue Theta condition turns the end face by the spatial
    // rotation vector Lambdaf & (Theta - Theta.prevIter())
    Theta_.boundaryFieldRef()[patchI][0] =
        Theta_.prevIter().boundaryField()[patchI][0]
      + (Lambdaf_.boundaryField()[patchI][0].T() & increment);
}


vector coupledTotalLagNewtonRaphsonBeam::rigidBodyEndBeamMoment() const
{
    if (!rigidBodyEndPtr_().welded())
    {
        return vector::zero;
    }

    // As for the force, the beam acts on the body with minus the end-face
    // moment
    return -M_.boundaryField()[rigidBodyEndPatchIndex_][0];
}

void coupledTotalLagNewtonRaphsonBeam::checkRigidBodyEndPinMoment() const
{
    const label patchI = rigidBodyEndPatchIndex_;

    // The moment the condition applied in the latest beam solve (a
    // momentSeries is evaluated when the condition is updated)
    const momentBeamRotationNRFvPatchVectorField& ThetaPatch =
        refCast<const momentBeamRotationNRFvPatchVectorField>
        (
            Theta_.boundaryField()[patchI]
        );

    if (max(mag(ThetaPatch.moment())) > SMALL)
    {
        FatalErrorInFunction
            << "A pinned rigidBodyEnd passes no moment: the "
            << "momentBeamRotationNR Theta condition on patch "
            << mesh().boundary()[patchI].name() << " must have zero moment "
            << "(also from a momentSeries), got " << ThetaPatch.moment()
            << " at time " << runTime().timeName() << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
