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

\*---------------------------------------------------------------------------*/

#include "rigidBodyEnd.H"
#include "Pstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rigidBodyEnd, 0);
}


// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * //

namespace
{
    // Skew tensor with spin(a) & b = a ^ b
    Foam::tensor spin(const Foam::vector& a)
    {
        return Foam::tensor
        (
            0,       -a.z(),  a.y(),
            a.z(),    0,     -a.x(),
           -a.y(),    a.x(),  0
        );
    }

    // Rotation tensor for a rotation vector (Rodrigues)
    Foam::tensor expMap(const Foam::vector& theta)
    {
        const Foam::scalar angle = Foam::mag(theta);
        const Foam::tensor K = spin(theta);

        Foam::scalar a, b;
        if (angle < 1e-6)
        {
            a = 1 - Foam::sqr(angle)/6;
            b = 0.5 - Foam::sqr(angle)/24;
        }
        else
        {
            a = Foam::sin(angle)/angle;
            b = (1 - Foam::cos(angle))/Foam::sqr(angle);
        }

        return Foam::tensor::I + a*K + b*(K & K);
    }

    // Tensor with a diagonal
    Foam::tensor diagonal(const Foam::diagTensor& d)
    {
        return Foam::tensor(d.xx(), 0, 0, 0, d.yy(), 0, 0, 0, d.zz());
    }

    // Rotation vector of a rotation tensor
    Foam::vector logMap(const Foam::tensor& R)
    {
        const Foam::vector skew
        (
            R.zy() - R.yz(),
            R.xz() - R.zx(),
            R.yx() - R.xy()
        );

        const Foam::scalar c =
            Foam::min(Foam::max(0.5*(Foam::tr(R) - 1), -1.0), 1.0);
        const Foam::scalar angle = Foam::acos(c);

        if (angle < 1e-6)
        {
            return 0.5*skew;
        }

        return angle/(2*Foam::sin(angle))*skew;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rigidBodyEnd::makeHistoryFile()
{
    if (!Pstream::master())
    {
        return;
    }

    const fileName historyDir
    (
        runTime_.path()/"postProcessing"/"rigidBodyEnd"/runTime_.timeName()
    );

    mkDir(historyDir);

    historyFilePtr_.reset(new OFstream(historyDir/"rigidBodyEnd.dat"));

    historyFilePtr_()
        << "# Time"
        << " dispX dispY dispZ"
        << " velX velY velZ"
        << " accX accY accZ"
        << " beamForceX beamForceY beamForceZ"
        << " couplingIterations"
        << " attachDispX attachDispY attachDispZ"
        << " rotationX rotationY rotationZ"
        << " omegaX omegaY omegaZ"
        << " beamTorqueX beamTorqueY beamTorqueZ"
        << endl;
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


Foam::vector Foam::rigidBodyEnd::velocityAt(const vector& x) const
{
    const scalar deltaT = runTime_.deltaTValue();
    const vector a = (x - predictor())/(beta_*sqr(deltaT));

    return v0_ + deltaT*((1 - gamma_)*a0_ + gamma_*a);
}


Foam::tensor Foam::rigidBodyEnd::inertia(const vector& theta) const
{
    const tensor R = orientation(theta);

    return (R & diagonal(momentOfInertia_) & R.T());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rigidBodyEnd::rigidBodyEnd
(
    const Time& runTime,
    const dictionary& dict,
    const word& defaultPatchName,
    const vector& gravity
)
:
    runTime_(runTime),
    patchName_(dict.getOrDefault<word>("patch", defaultPatchName)),
    mass_(dict.get<scalar>("mass")),
    momentOfInertia_
    (
        dict.getOrDefault<diagTensor>("momentOfInertia", diagTensor::zero)
    ),
    rotation_(cmptMax(momentOfInertia_) > SMALL),
    centreOfMassGiven_(dict.found("centreOfMass")),
    centreOfMass_(dict.getOrDefault<point>("centreOfMass", point::zero)),
    arm0_(vector::zero),
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
    nCouplingIterations_(dict.getOrDefault<label>("nCouplingIterations", 1)),
    couplingTolerance_(dict.getOrDefault<scalar>("couplingTolerance", 1e-8)),
    couplingRelaxation_(dict.getOrDefault<scalar>("couplingRelaxation", 1)),
    nJacobianChecks_(dict.getOrDefault<label>("jacobianCheck", 0)),
    x_(vector::zero),
    v_(dict.getOrDefault<vector>("velocity", vector::zero)),
    a_(vector::zero),
    x0_(vector::zero),
    v0_(vector::zero),
    a0_(vector::zero),
    Q_(tensor::I),
    omega_(dict.getOrDefault<vector>("angularVelocity", vector::zero)),
    alpha_(vector::zero),
    Q0_(tensor::I),
    omega0_(vector::zero),
    alpha0_(vector::zero),
    theta_(vector::zero),
    beamForce_(vector::zero),
    beamTorque_(vector::zero),
    attachmentOffset_(vector::zero),
    timeIndex_(runTime.timeIndex()),
    nIterationsUsed_(0),
    historyFilePtr_()
{
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

    if (!rotation_)
    {
        omega_ = vector::zero;
    }

    // The beam starts unloaded, so the initial accelerations come from the
    // applied loads, gravity and damping only
    a_ = (totalForce(vector::zero) - linearDamping_*v_)/mass_;

    if (rotation_)
    {
        const tensor J = diagonal(momentOfInertia_);
        alpha_ =
            inv(J)
          & (
                externalMoment_
              - (omega_ ^ (J & omega_))
              - angularDamping_*omega_
            );
    }

    Info<< "rigidBodyEnd on patch " << patchName_ << nl
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

    makeHistoryFile();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rigidBodyEnd::setAttachmentPoint(const point& attachmentPoint)
{
    if (!centreOfMassGiven_)
    {
        centreOfMass_ = attachmentPoint;
    }

    arm0_ = attachmentPoint - centreOfMass_;

    if (!rotation_ && mag(arm0_) > SMALL)
    {
        FatalErrorInFunction
            << "rigidBodyEnd centreOfMass " << centreOfMass_
            << " is away from the attachment point " << attachmentPoint
            << ": give a momentOfInertia so the body can rotate"
            << abort(FatalError);
    }

    Info<< "rigidBodyEnd: attachment point " << attachmentPoint
        << ", centre of mass " << centreOfMass_
        << ", arm " << arm0_ << endl;
}


Foam::vector Foam::rigidBodyEnd::totalForce(const vector& beamForce) const
{
    return beamForce + externalForce_ + mass_*gravity_;
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


Foam::vector Foam::rigidBodyEnd::predictor() const
{
    const scalar deltaT = runTime_.deltaTValue();

    return x0_ + deltaT*v0_ + sqr(deltaT)*(0.5 - beta_)*a0_;
}


Foam::scalar Foam::rigidBodyEnd::inertiaCoefficient() const
{
    return mass_/(beta_*sqr(runTime_.deltaTValue()));
}


Foam::vector Foam::rigidBodyEnd::translationalResidual
(
    const vector& x,
    const vector& beamForce
) const
{
    return
        inertiaCoefficient()*(x - predictor())
      + linearDamping_*velocityAt(x)
      - totalForce(beamForce);
}


Foam::scalar Foam::rigidBodyEnd::translationalJacobian() const
{
    const scalar deltaT = runTime_.deltaTValue();

    return inertiaCoefficient() + linearDamping_*gamma_/(beta_*deltaT);
}


Foam::tensor Foam::rigidBodyEnd::orientation(const vector& theta) const
{
    return (expMap(theta) & Q0_);
}


Foam::vector Foam::rigidBodyEnd::arm(const vector& theta) const
{
    return (orientation(theta) & arm0_);
}


Foam::vector Foam::rigidBodyEnd::attachmentDisplacement
(
    const vector& x,
    const vector& theta
) const
{
    return x + arm(theta) - arm0_;
}


Foam::vector Foam::rigidBodyEnd::rotationalResidual
(
    const vector& theta,
    const vector& beamForce
) const
{
    if (!rotation_)
    {
        return vector::zero;
    }

    const tensor J = inertia(theta);
    const vector omega = angularVelocity(theta);

    return
        (J & angularAcceleration(theta))
      + (omega ^ (J & omega))
      + angularDamping_*omega
      - externalMoment_
      - (arm(theta) ^ beamForce);
}


Foam::tensor Foam::rigidBodyEnd::rotationalInertiaJacobian
(
    const vector& theta
) const
{
    if (!rotation_)
    {
        return tensor::zero;
    }

    const scalar deltaT = runTime_.deltaTValue();
    const tensor J = inertia(theta);
    const vector omega = angularVelocity(theta);
    const vector alpha = angularAcceleration(theta);

    // d(J alpha): alpha changes by dtheta/(beta dt^2), and J turns with the
    // body by rotationJacobian & dtheta; d(omega x J omega) and d(cr omega):
    // omega changes by gamma/(beta dt) dtheta
    return
        J/(beta_*sqr(deltaT))
      + ((-spin(J & alpha) + (J & spin(alpha))) & rotationJacobian(theta))
      + gamma_/(beta_*deltaT)
       *((spin(omega) & J) - spin(J & omega) + angularDamping_*tensor::I);
}


Foam::tensor Foam::rigidBodyEnd::rotationJacobian(const vector& theta) const
{
    const scalar angle = mag(theta);
    const tensor K = spin(theta);

    scalar a, b;
    if (angle < 1e-4)
    {
        a = 0.5 - sqr(angle)/24;
        b = 1.0/6.0 - sqr(angle)/120;
    }
    else
    {
        a = (1 - cos(angle))/sqr(angle);
        b = (angle - sin(angle))/pow3(angle);
    }

    return tensor::I + a*K + b*(K & K);
}


Foam::vector Foam::rigidBodyEnd::newmarkRotation
(
    const vector& beamForce,
    const vector& thetaEstimate
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
          - (omega ^ (J & omega))
          - angularDamping_*(omega0_ + deltaT*(1 - gamma_)*alpha0_)
        );

    return
        angularPredictor()
      + beta_*sqr(runTime_.deltaTValue())*alphaNew;
}


void Foam::rigidBodyEnd::accept
(
    const vector& displacement,
    const vector& theta,
    const vector& beamForce,
    const label nIterations
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
    }

    beamForce_ = beamForce;
    beamTorque_ = (arm(theta_) ^ beamForce);
    nIterationsUsed_ = nIterations;
}


void Foam::rigidBodyEnd::write() const
{
    if (!historyFilePtr_.valid())
    {
        return;
    }

    OFstream& os = historyFilePtr_();

    const vector attach = attachmentDisplacement(x_, theta_);
    const vector rotation = logMap(Q_);

    os  << runTime_.value()
        << " " << x_.x() << " " << x_.y() << " " << x_.z()
        << " " << v_.x() << " " << v_.y() << " " << v_.z()
        << " " << a_.x() << " " << a_.y() << " " << a_.z()
        << " " << beamForce_.x()
        << " " << beamForce_.y()
        << " " << beamForce_.z()
        << " " << nIterationsUsed_
        << " " << attach.x() << " " << attach.y() << " " << attach.z()
        << " " << rotation.x() << " " << rotation.y() << " " << rotation.z()
        << " " << omega_.x() << " " << omega_.y() << " " << omega_.z()
        << " " << beamTorque_.x()
        << " " << beamTorque_.y()
        << " " << beamTorque_.z()
        << endl;
}


// ************************************************************************* //
