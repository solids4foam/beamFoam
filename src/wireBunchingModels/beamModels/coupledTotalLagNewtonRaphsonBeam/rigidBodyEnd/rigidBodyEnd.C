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

\*---------------------------------------------------------------------------*/

#include "rigidBodyEnd.H"
#include "pseudoVector.H"
#include "Pstream.H"
#include "OSspecific.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rigidBodyEnd, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rigidBodyEnd::readState()
{
    x_ = state_.get<vector>("displacement");
    v_ = state_.get<vector>("velocity");
    a_ = state_.get<vector>("acceleration");
    Q_ = state_.get<tensor>("orientation");
    omega_ = state_.get<vector>("angularVelocity");
    alpha_ = state_.get<vector>("angularAcceleration");
    Qinit_ = state_.get<tensor>("initialOrientation");
    rotationSum_ = state_.get<vector>("rotationSum");
    beamForce_ = state_.get<vector>("beamForce");
    beamMoment_ = state_.get<vector>("beamMoment");
    beamTorque_ = state_.get<vector>("beamTorque");
    attachmentOffset_ = state_.get<vector>("attachmentOffset");
    nNotConverged_ = state_.get<label>("nNotConverged");
}


void Foam::rigidBodyEnd::storeState()
{
    state_.set("displacement", x_);
    state_.set("velocity", v_);
    state_.set("acceleration", a_);
    state_.set("orientation", Q_);
    state_.set("angularVelocity", omega_);
    state_.set("angularAcceleration", alpha_);
    state_.set("initialOrientation", Qinit_);
    state_.set("rotationSum", rotationSum_);
    state_.set("beamForce", beamForce_);
    state_.set("beamMoment", beamMoment_);
    state_.set("beamTorque", beamTorque_);
    state_.set("attachmentOffset", attachmentOffset_);
    state_.set("nNotConverged", nNotConverged_);
}


void Foam::rigidBodyEnd::makeHistoryFile(const objectRegistry& db)
{
    if (!Pstream::master())
    {
        return;
    }

    // Named after the start time, so a restarted run writes a new file. As
    // for function objects, a beam that is not in the default region writes
    // under postProcessing/<region>, so bodies in different regions do not
    // share a file
    fileName historyDir
    (
        runTime_.path()/"postProcessing"/polyMesh::regionName(db.name())
       /"rigidBodyEnd"/runTime_.timeName(runTime_.startTime().value())
    );
    historyDir.clean();

    mkDir(historyDir);

    historyFilePtr_.reset(new OFstream(historyDir/"rigidBodyEnd.dat"));

    historyFilePtr_()
        << "# Time"
        << " dispX dispY dispZ"
        << " velX velY velZ"
        << " accX accY accZ"
        << " beamForceX beamForceY beamForceZ"
        << " couplingIterations couplingConverged nNotConverged"
        << " attachDispX attachDispY attachDispZ"
        << " Qxx Qxy Qxz Qyx Qyy Qyz Qzx Qzy Qzz"
        << " rotationSumX rotationSumY rotationSumZ"
        << " omegaX omegaY omegaZ"
        << " beamTorqueX beamTorqueY beamTorqueZ"
        << " beamMomentX beamMomentY beamMomentZ"
        << endl;
}


Foam::vector Foam::rigidBodyEnd::predictor() const
{
    const scalar deltaT = runTime_.deltaTValue();

    return x0_ + deltaT*v0_ + sqr(deltaT)*(0.5 - beta_)*a0_;
}


Foam::vector Foam::rigidBodyEnd::angularPredictor() const
{
    const scalar deltaT = runTime_.deltaTValue();

    return deltaT*omega0_ + sqr(deltaT)*(0.5 - beta_)*alpha0_;
}


Foam::vector Foam::rigidBodyEnd::angularAcceleration
(
    const vector& theta
) const
{
    return (theta - angularPredictor())/(beta_*sqr(runTime_.deltaTValue()));
}


Foam::vector Foam::rigidBodyEnd::angularVelocity(const vector& theta) const
{
    const scalar deltaT = runTime_.deltaTValue();

    return
        omega0_
      + deltaT*((1 - gamma_)*alpha0_ + gamma_*angularAcceleration(theta));
}


Foam::tensor Foam::rigidBodyEnd::inertia(const vector& theta) const
{
    const tensor R = orientation(theta);

    return (R & momentOfInertia_ & R.T());
}


Foam::tensor Foam::rigidBodyEnd::orientation(const vector& theta) const
{
    return (stepRotation(theta) & Q0_);
}


Foam::vector Foam::rigidBodyEnd::totalForce(const vector& beamForce) const
{
    return beamForce + externalForce_ + mass_*gravity_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyEnd::rigidBodyEnd
(
    const objectRegistry& db,
    const dictionary& dict,
    const word& defaultPatchName,
    const vector& gravity
)
:
    runTime_(db.time()),
    patchName_(dict.getOrDefault<word>("patch", defaultPatchName)),
    joint_(jointType::pin),
    mass_(dict.get<scalar>("mass")),
    momentOfInertia_
    (
        dict.getOrDefault<diagTensor>("momentOfInertia", diagTensor::zero)
    ),
    rotation_(cmptMax(momentOfInertia_) > SMALL),
    centreOfMassGiven_(dict.found("centreOfMass")),
    centreOfMass_(dict.getOrDefault<point>("centreOfMass", point::zero)),
    arm0_(vector::zero),
    armBody_(vector::zero),
    externalForce_(dict.getOrDefault<vector>("externalForce", vector::zero)),
    externalMoment_
    (
        dict.getOrDefault<vector>("externalMoment", vector::zero)
    ),
    linearDamping_(dict.getOrDefault<scalar>("linearDamping", 0)),
    angularDamping_(dict.getOrDefault<scalar>("angularDamping", 0)),
    gravity_(gravity),
    beta_(dict.getOrDefault<scalar>("newmarkBeta", 0.25)),
    gamma_(dict.getOrDefault<scalar>("newmarkGamma", 0.5)),
    nCouplingIterations_
    (
        dict.getOrDefault<label>("nCouplingIterations", 20)
    ),
    couplingTolerance_(dict.getOrDefault<scalar>("couplingTolerance", 1e-8)),
    couplingRelaxation_(dict.getOrDefault<scalar>("couplingRelaxation", 1)),
    state_
    (
        IOobject
        (
            "rigidBodyEnd",
            runTime_.timeName(),
            "uniform",
            db,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    restarted_(state_.found("displacement")),
    x_(vector::zero),
    v_(dict.getOrDefault<vector>("velocity", vector::zero)),
    a_(vector::zero),
    x0_(vector::zero),
    v0_(vector::zero),
    a0_(vector::zero),
    Q_(dict.getOrDefault<tensor>("orientation", tensor::I)),
    omega_(dict.getOrDefault<vector>("angularVelocity", vector::zero)),
    alpha_(vector::zero),
    Qinit_(Q_),
    Q0_(Q_),
    omega0_(vector::zero),
    alpha0_(vector::zero),
    theta_(vector::zero),
    rotationSum_(vector::zero),
    rotationSum0_(vector::zero),
    beamForce_(vector::zero),
    beamMoment_(vector::zero),
    beamTorque_(vector::zero),
    attachmentOffset_(vector::zero),
    timeIndex_(runTime_.timeIndex()),
    nIterationsUsed_(0),
    converged_(true),
    nNotConverged_(0),
    historyFilePtr_()
{
    const word joint(dict.getOrDefault<word>("joint", "pin"));

    if (joint == "welded")
    {
        joint_ = jointType::welded;
    }
    else if (joint != "pin")
    {
        FatalIOErrorInFunction(dict)
            << "Invalid rigidBodyEnd joint " << joint << nl
            << "Valid options are pin and welded"
            << exit(FatalIOError);
    }

    if (mass_ <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "rigidBodyEnd mass must be positive, got " << mass_
            << exit(FatalIOError);
    }

    if (nCouplingIterations_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "nCouplingIterations must be at least 1"
            << exit(FatalIOError);
    }

    if (rotation_ && cmptMin(momentOfInertia_) <= SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "All principal moments of inertia must be positive when "
            << "rotation is modelled, got " << momentOfInertia_
            << exit(FatalIOError);
    }

    if (welded() && !rotation_)
    {
        FatalIOErrorInFunction(dict)
            << "A welded rigidBodyEnd turns the beam end with the body: "
            << "give a momentOfInertia" << exit(FatalIOError);
    }

    if (!rotation_)
    {
        omega_ = vector::zero;
    }

    if (restarted_)
    {
        readState();

        Info<< "rigidBodyEnd: state read from " << state_.objectPath()
            << endl;
    }

    Info<< "rigidBodyEnd on patch " << patchName_ << nl
        << "    joint " << joint << nl
        << "    mass " << mass_ << nl
        << "    momentOfInertia " << momentOfInertia_
        << (rotation_ ? "" : " (point mass: rotation not modelled)") << nl
        << "    externalForce " << externalForce_
        << ", externalMoment " << externalMoment_ << nl
        << "    linearDamping " << linearDamping_
        << ", angularDamping " << angularDamping_ << nl
        << "    gravity " << gravity_ << nl
        << "    Newmark beta " << beta_ << ", gamma " << gamma_ << nl
        << "    nCouplingIterations " << nCouplingIterations_
        << ", couplingTolerance " << couplingTolerance_
        << ", couplingRelaxation " << couplingRelaxation_ << endl;

    makeHistoryFile(db);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rigidBodyEnd::attach
(
    const point& endPoint,
    const vector& endW,
    const vector& beamForce,
    const vector& beamMoment
)
{
    // At the start of a run the body is at zero displacement, so the current
    // end displacement (e.g. the placement from setInitialPositionBeam or a
    // beamInit state) is the offset. On a restart W includes the body motion
    // and the offset comes from the saved state
    if (!restarted_)
    {
        attachmentOffset_ = endW;
    }

    const point attachmentPoint(endPoint + attachmentOffset_);

    if (!centreOfMassGiven_)
    {
        centreOfMass_ = attachmentPoint;
    }

    arm0_ = attachmentPoint - centreOfMass_;
    armBody_ = (Qinit_.T() & arm0_);

    if (!rotation_ && mag(arm0_) > SMALL)
    {
        FatalErrorInFunction
            << "rigidBodyEnd centreOfMass " << centreOfMass_
            << " is away from the attachment point " << attachmentPoint
            << ": give a momentOfInertia so the body can rotate"
            << abort(FatalError);
    }

    // The old-time state is the current state, so that the first time step
    // starts from it
    Q0_ = Q_;
    rotationSum0_ = rotationSum_;

    if (!restarted_)
    {
        // The beam may already be loaded (a pretensioned line, a case
        // started from a beamInit state), so the initial accelerations
        // include the current force and moment from the beam
        beamForce_ = beamForce;
        beamMoment_ = welded() ? beamMoment : vector::zero;
        beamTorque_ = (arm(vector::zero) ^ beamForce_) + beamMoment_;

        a_ = (totalForce(beamForce_) - linearDamping_*v_)/mass_;

        if (rotation_)
        {
            const tensor J = inertia(vector::zero);

            alpha_ =
                inv(J)
              & (
                    externalMoment_
                  + beamTorque_
                  - (omega_ ^ (J & omega_))
                  - angularDamping_*omega_
                );
        }
    }

    x0_ = x_;
    v0_ = v_;
    a0_ = a_;
    omega0_ = omega_;
    alpha0_ = alpha_;

    storeState();

    Info<< "rigidBodyEnd: attachment point " << attachmentPoint
        << ", centre of mass " << centreOfMass_
        << ", arm " << arm0_ << nl
        << "    initial beam force " << beamForce_
        << ", beam moment " << beamMoment_ << nl
        << "    acceleration " << a_
        << ", angular acceleration " << alpha_ << endl;
}


void Foam::rigidBodyEnd::newTimeStep()
{
    if (runTime_.timeIndex() == timeIndex_)
    {
        return;
    }

    timeIndex_ = runTime_.timeIndex();

    x0_ = x_;
    v0_ = v_;
    a0_ = a_;

    Q0_ = Q_;
    omega0_ = omega_;
    alpha0_ = alpha_;
    theta_ = vector::zero;
    rotationSum0_ = rotationSum_;

    // Not yet counted as a time step that did not converge
    converged_ = true;
}


Foam::vector Foam::rigidBodyEnd::newmarkDisplacement
(
    const vector& beamForce
) const
{
    const scalar deltaT = runTime_.deltaTValue();

    // m a + c (v0 + dt((1 - gamma) a0 + gamma a)) = total force
    const vector aNew =
        (
            totalForce(beamForce)
          - linearDamping_*(v0_ + deltaT*(1 - gamma_)*a0_)
        )
       /(mass_ + linearDamping_*gamma_*deltaT);

    return
        x0_
      + deltaT*v0_
      + sqr(deltaT)*((0.5 - beta_)*a0_ + beta_*aNew);
}


Foam::vector Foam::rigidBodyEnd::newmarkRotation
(
    const vector& beamForce,
    const vector& thetaEstimate,
    const vector& beamMoment
) const
{
    if (!rotation_)
    {
        return vector::zero;
    }

    const scalar deltaT = runTime_.deltaTValue();
    const tensor J = inertia(thetaEstimate);
    const vector omega = angularVelocity(thetaEstimate);

    // Damping taken implicitly, the gyroscopic term from the estimate
    const vector alphaNew =
        inv(J + angularDamping_*gamma_*deltaT*tensor::I)
      & (
            externalMoment_
          + (arm(thetaEstimate) ^ beamForce)
          + beamMoment
          - (omega ^ (J & omega))
          - angularDamping_*(omega0_ + deltaT*(1 - gamma_)*alpha0_)
        );

    return angularPredictor() + beta_*sqr(deltaT)*alphaNew;
}


Foam::scalar Foam::rigidBodyEnd::totalRotation(const vector& theta) const
{
    return mag(rotationSum0_ + theta);
}


Foam::vector Foam::rigidBodyEnd::arm(const vector& theta) const
{
    return (orientation(theta) & armBody_);
}


Foam::vector Foam::rigidBodyEnd::attachmentDisplacement
(
    const vector& x,
    const vector& theta
) const
{
    return x + arm(theta) - arm0_;
}


Foam::tensor Foam::rigidBodyEnd::stepRotation(const vector& theta) const
{
    return rotationMatrix(theta);
}


void Foam::rigidBodyEnd::accept
(
    const vector& displacement,
    const vector& theta,
    const vector& beamForce,
    const vector& beamMoment,
    const label nIterations,
    const bool converged
)
{
    const scalar deltaT = runTime_.deltaTValue();

    x_ = displacement;
    a_ = (x_ - predictor())/(beta_*sqr(deltaT));
    v_ = v0_ + deltaT*((1 - gamma_)*a0_ + gamma_*a_);

    if (rotation_)
    {
        theta_ = theta;
        Q_ = orientation(theta);
        alpha_ = angularAcceleration(theta);
        omega_ = angularVelocity(theta);
        rotationSum_ = rotationSum0_ + theta;
    }

    beamForce_ = beamForce;
    beamMoment_ = beamMoment;
    beamTorque_ = (arm(theta_) ^ beamForce) + beamMoment;
    nIterationsUsed_ = nIterations;

    if (!converged && converged_)
    {
        ++nNotConverged_;
    }
    else if (converged && !converged_)
    {
        // The step was accepted before (another evolve() in the same time
        // step) without converging: it has now converged
        --nNotConverged_;
    }
    converged_ = converged;

    storeState();
}


void Foam::rigidBodyEnd::write() const
{
    if (!historyFilePtr_)
    {
        return;
    }

    OFstream& os = historyFilePtr_();

    const vector attach = attachmentDisplacement(x_, theta_);

    os  << runTime_.value()
        << " " << x_.x() << " " << x_.y() << " " << x_.z()
        << " " << v_.x() << " " << v_.y() << " " << v_.z()
        << " " << a_.x() << " " << a_.y() << " " << a_.z()
        << " " << beamForce_.x()
        << " " << beamForce_.y()
        << " " << beamForce_.z()
        << " " << nIterationsUsed_
        << " " << converged_
        << " " << nNotConverged_
        << " " << attach.x() << " " << attach.y() << " " << attach.z()
        << " " << Q_.xx() << " " << Q_.xy() << " " << Q_.xz()
        << " " << Q_.yx() << " " << Q_.yy() << " " << Q_.yz()
        << " " << Q_.zx() << " " << Q_.zy() << " " << Q_.zz()
        << " " << rotationSum_.x()
        << " " << rotationSum_.y()
        << " " << rotationSum_.z()
        << " " << omega_.x() << " " << omega_.y() << " " << omega_.z()
        << " " << beamTorque_.x()
        << " " << beamTorque_.y()
        << " " << beamTorque_.z()
        << " " << beamMoment_.x()
        << " " << beamMoment_.y()
        << " " << beamMoment_.z()
        << endl;
}


// ************************************************************************* //
