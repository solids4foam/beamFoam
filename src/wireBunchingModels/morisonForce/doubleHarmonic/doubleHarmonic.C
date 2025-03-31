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

#include "doubleHarmonic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(doubleHarmonic, 0);
addToRunTimeSelectionTable(morisonForce, doubleHarmonic, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


doubleHarmonic::doubleHarmonic
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& sd
)
:
    morisonForce(runTime, mesh, sd),

    Um1_
    (
        this->subDict(typeName + "Coeffs").lookup("Um1")
    ),

    Um2_
    (
        this->subDict(typeName + "Coeffs").lookup("Um2")
    ),

    k_
    (
        this->subDict(typeName + "Coeffs").lookup("waveNumber")
    ),

    period_
    (
        readScalar(this->subDict(typeName + "Coeffs").lookup("period"))
    ),

    phi_
    (
        readScalar(this->subDict(typeName + "Coeffs").lookup("phi"))
    ),

    omega_
    (
        2.0*M_PI/period_
    ),

    Tsoft_
    (
        this->subDict(typeName + "Coeffs").lookupOrDefault<scalar>("Tsoft", period_)
    )
{
}


void doubleHarmonic::printCoeffs()
{
    Info << "Loading Morison force: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector doubleHarmonic::U
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    return dimensionedVector
        (
            "null", dimVelocity,
            factor*Um1_*Foam::cos(omega_*t - (k_ & x) + phi_)
          + factor*Um2_*Foam::cos(2.0*(omega_*t - (k_ & x) + phi_))
        );
}


dimensionedVector doubleHarmonic::acc
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    return dimensionedVector
        (
            "null", dimVelocity/dimTime,
          - factor*omega_*Um1_*Foam::sin(omega_*t - (k_ & x) + phi_)
          - factor*2*omega_*Um2_*Foam::sin(2.0*(omega_*t - (k_ & x) + phi_))
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
