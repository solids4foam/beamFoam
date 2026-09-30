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
        << endl;
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
    externalForce_(dict.getOrDefault<vector>("externalForce", vector::zero)),
    gravity_(gravity),
    beta_(dict.getOrDefault<scalar>("newmarkBeta", 0.25)),
    gamma_(dict.getOrDefault<scalar>("newmarkGamma", 0.5)),
    nCouplingIterations_(dict.getOrDefault<label>("nCouplingIterations", 1)),
    couplingTolerance_(dict.getOrDefault<scalar>("couplingTolerance", 1e-8)),
    couplingRelaxation_(dict.getOrDefault<scalar>("couplingRelaxation", 1)),
    nJacobianChecks_(dict.getOrDefault<label>("jacobianCheck", 0)),
    x_(vector::zero),
    v_(vector::zero),
    a_(vector::zero),
    x0_(vector::zero),
    v0_(vector::zero),
    a0_(vector::zero),
    beamForce_(vector::zero),
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

    // The beam starts unloaded, so the initial acceleration comes from the
    // applied force and gravity only
    a_ = totalForce(vector::zero)/mass_;

    Info<< "rigidBodyEnd on patch " << patchName_ << nl
        << "    mass " << mass_ << nl
        << "    externalForce " << externalForce_ << nl
        << "    gravity " << gravity_ << nl
        << "    Newmark beta " << beta_ << ", gamma " << gamma_ << nl
        << "    nCouplingIterations " << nCouplingIterations_
        << ", couplingTolerance " << couplingTolerance_
        << ", couplingRelaxation " << couplingRelaxation_ << nl
        << "    point mass: rotation not modelled" << endl;

    makeHistoryFile();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
}


Foam::vector Foam::rigidBodyEnd::newmarkDisplacement
(
    const vector& beamForce
) const
{
    const scalar deltaT = runTime_.deltaTValue();
    const vector aNew = totalForce(beamForce)/mass_;

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


void Foam::rigidBodyEnd::accept
(
    const vector& displacement,
    const vector& beamForce,
    const label nIterations
)
{
    const scalar deltaT = runTime_.deltaTValue();

    x_ = displacement;
    a_ = (x_ - predictor())/(beta_*sqr(deltaT));
    v_ = v0_ + deltaT*((1 - gamma_)*a0_ + gamma_*a_);

    beamForce_ = beamForce;
    nIterationsUsed_ = nIterations;
}


void Foam::rigidBodyEnd::write() const
{
    if (!historyFilePtr_.valid())
    {
        return;
    }

    OFstream& os = historyFilePtr_();

    os  << runTime_.value()
        << " " << x_.x() << " " << x_.y() << " " << x_.z()
        << " " << v_.x() << " " << v_.y() << " " << v_.z()
        << " " << a_.x() << " " << a_.y() << " " << a_.z()
        << " " << beamForce_.x()
        << " " << beamForce_.y()
        << " " << beamForce_.z()
        << " " << nIterationsUsed_
        << endl;
}


// ************************************************************************* //
