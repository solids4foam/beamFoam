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

#include "streamFunctionUserCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(streamFunctionUserCoeffs, 0);
addToRunTimeSelectionTable(morisonForce, streamFunctionUserCoeffs, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


streamFunctionUserCoeffs::streamFunctionUserCoeffs
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& sd
)
:
    morisonForce(runTime, mesh, sd),

    N_(readLabel(this->subDict(typeName + "Coeffs").lookup("N"))),

    omega_("omega", this->subDict(typeName + "Coeffs"), N_),
    velAmp_("velAmp", this->subDict(typeName + "Coeffs"), N_),
    waveNumber_("waveNumber", this->subDict(typeName + "Coeffs"), N_),
    phases_("phases", this->subDict(typeName + "Coeffs"), N_),

    meanVel_(readScalar(this->subDict(typeName + "Coeffs").lookup("meanVel"))),

    h_(readScalar(this->subDict(typeName + "Coeffs").lookup("depth"))),

    Tsoft_(readScalar(this->subDict(typeName + "Coeffs").lookup("Tsoft"))),

    direction_(this->subDict(typeName + "Coeffs").lookup("direction"))
{
}

void streamFunctionUserCoeffs::printCoeffs()
{
    Info << "Loading Morison force: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector streamFunctionUserCoeffs::U
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    // Direct elimination of water depth
    scalar Z(x & vertDir_);
    scalar Uhorz(meanVel_);
    scalar Uvert(0);

    const scalarField arg(omega_*t + phases_ - waveNumber_*(direction_ & x));

    forAll (arg, cmpi)
    {
        Uhorz += velAmp_[cmpi]*Foam::cosh(waveNumber_[cmpi]*Z)*Foam::cos(arg[cmpi]);
        Uvert -= velAmp_[cmpi]*Foam::sinh(waveNumber_[cmpi]*Z)*Foam::sin(arg[cmpi]);
    }

    Uhorz *= factor;
    Uvert *= factor;

    vector vel(Uhorz*direction_);

    vel += Uvert*vertDir_;

    return dimensionedVector("null", dimVelocity, vel);
}


dimensionedVector streamFunctionUserCoeffs::acc
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    // Direct elimination of water depth
    scalar Z(x & vertDir_);
    scalar Ahorz(0);
    scalar Avert(0);

    const scalarField arg(omega_*t + phases_ - waveNumber_*(direction_ & x));

    forAll (arg, cmpi)
    {
        Ahorz += -omega_[cmpi]*velAmp_[cmpi]*Foam::cosh(waveNumber_[cmpi]*Z)*Foam::sin(arg[cmpi]);
        Avert += -omega_[cmpi]*velAmp_[cmpi]*Foam::sinh(waveNumber_[cmpi]*Z)*Foam::cos(arg[cmpi]);
    }

    Ahorz *= factor;
    Avert *= factor;

    vector acc(Ahorz*direction_);

    acc += Avert*vertDir_;

    return dimensionedVector("null", dimVelocity/dimTime, acc);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
