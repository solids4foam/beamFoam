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

#include "userDefinedWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(userDefinedWave, 0);
addToRunTimeSelectionTable(morisonForce, userDefinedWave, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


userDefinedWave::userDefinedWave
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& sd
)
:
    morisonForce(runTime, mesh, sd),

    includeVertical_(Switch(this->subDict(typeName + "Coeffs").lookup("includeVertical"))),
    includeWaveNumber_(Switch(this->subDict(typeName + "Coeffs").lookup("includeWaveNumber")))
{
    setCoefficients();

    spaceFact_ = 0.0;
    if (includeWaveNumber_)
    {
        spaceFact_ = 1.0;
    }

}


void userDefinedWave::setCoefficients()
{
    const dictionary sd(this->subDict(typeName + "Coeffs"));

    N_ = readLabel(sd.lookup("N"));

    // Read the velocity coefficients
    u0_ = readScalar(sd.lookup("U0"));
    du0dz_ = readScalar(sd.lookup("dU0dz"));


    {
        scalarField coeffs("uSin", sd, N_);
        uSin_.setSize(coeffs.size());
        uSin_ = coeffs;
    }

    {
        scalarField coeffs("uCos", sd, N_);
        uCos_.setSize(coeffs.size());
        uCos_ = coeffs;
    }

    {
        scalarField coeffs("wSin", sd, N_);
        wSin_.setSize(coeffs.size());
        wSin_ = coeffs;
    }

    {
        scalarField coeffs("wCos", sd, N_);
        wCos_.setSize(coeffs.size());
        wCos_ = coeffs;
    }

    // Convert to the correct OpenFoam (dimensioned) format
    h_ = readScalar(sd.lookup("depth"));
    period_ = readScalar(sd.lookup("period"));

    omega_ = 2.0*M_PI/period_;

    Tsoft_ = period_;
    if (sd.found( "Tsoft" ))
    {
        Tsoft_ = readScalar(sd.lookup("Tsoft"));
    }

    vector direction( vector(sd.lookup("direction")));
    direction /= Foam::mag(direction);
    k_ = readScalar(sd.lookup("waveNumber"))*direction;
    K_ = Foam::mag(k_);
}


void userDefinedWave::printCoeffs()
{
    Info << "Loading Morison force: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector userDefinedWave::U
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    scalar Z((x & vertDir_) - h_);

    scalar Uhorz(u0_ + du0dz_*(x & vertDir_));
    scalar Uvert(0);

    scalar arg = omega_*t - (k_ & x)*spaceFact_;

    forAll (uSin_, ii)
    {
        scalar n = ii + 1;

        Uhorz += uSin_[ii]*Foam::cosh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::sin(n*arg);
        Uhorz += uCos_[ii]*Foam::cosh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::cos(n*arg);

        Uvert += wSin_[ii]*Foam::sinh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::sin(n*arg);
        Uvert += wCos_[ii]*Foam::sinh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::cos(n*arg);
    }


    Uhorz *= factor;
    Uvert *= factor;

    vector vel(Uhorz*k_/K_);

    if (includeVertical_)
    {
        vel += Uvert*vertDir_;
    }

    return dimensionedVector("null", dimVelocity, vel);
}


dimensionedVector userDefinedWave::acc
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    scalar Z((x & vertDir_) - h_);

    scalar ahorz(0);
    scalar avert(0);

    scalar arg = omega_*t - (k_ & x)*spaceFact_;

    forAll (uSin_, ii)
    {
        scalar n = ii + 1;

        ahorz += n*omega_*uSin_[ii]*Foam::cosh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::cos(n*arg);
        ahorz -= n*omega_*uCos_[ii]*Foam::cosh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::sin(n*arg);

        avert += n*omega_*wSin_[ii]*Foam::sinh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::cos(n*arg);
        avert -= n*omega_*wCos_[ii]*Foam::sinh(n*K_*(Z + h_))/Foam::cosh(n*K_*h_)
            *Foam::sin(n*arg);
    }

    ahorz *= factor;
    avert *= factor;

    vector acc(ahorz*k_/K_);

    if (includeVertical_)
    {
        acc += avert*vertDir_;
    }

    return dimensionedVector("null", dimVelocity/dimTime, acc);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
