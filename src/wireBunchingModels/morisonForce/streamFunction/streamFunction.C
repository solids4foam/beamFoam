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

#include "streamFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(streamFunction, 0);
addToRunTimeSelectionTable(morisonForce, streamFunction, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


streamFunction::streamFunction
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


void streamFunction::setCoefficients()
{
    writeInputFile();

    // Initialise the coefficients
    A_.setSize(N_, 0.0);
    B_.setSize(N_, 0.0);

    // Execute the fortran program 'fenton4Foam'. Name change not to collide
    // with other installations of the same code.
    system("fenton4Foam > /dev/null");



    // Read the output file
    IOdictionary fentonCoeffs
    (
        IOobject
        (
            "fenton4Foam",
            runTime_.constant(),
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read the various output variables
    dimensionedScalar depthNonDim(fentonCoeffs.lookup("depth"));
    dimensionedScalar periodNonDim(fentonCoeffs.lookup("period"));
    dimensionedScalar uBarNonDim(fentonCoeffs.lookup("uBar"));
    const dictionary sd(this->subDict(typeName + "Coeffs"));

    scalarField coeffsA("aCoeffs", fentonCoeffs, N_);
    scalarField coeffsB("bCoeffs", fentonCoeffs, N_);

    A_ = coeffsA;
    B_ = coeffsB;


    // Convert to the correct OpenFoam (dimensioned) format
    h_ = readScalar(sd.lookup("depth"));
    scalar waveNumber = depthNonDim.value()/h_;

    period_ = periodNonDim.value()/Foam::sqrt(G_*waveNumber);
    omega_ = 2.0*M_PI/period_;

    uBar_ = uBarNonDim.value()*Foam::sqrt(G_/waveNumber);
    A_ /= waveNumber;

    forAll(B_, coeffi)
    {
        B_[coeffi] *= (coeffi + 1)*Foam::sqrt(G_/waveNumber);
    }

    Tsoft_ = period_;
    if (sd.found( "Tsoft" ))
    {
        Tsoft_ = readScalar(sd.lookup("Tsoft"));
    }

    vector direction( vector(sd.lookup("direction")));
    direction /= Foam::mag(direction);
    k_ = waveNumber*direction;
}

void streamFunction::writeInputFile()
{
    const dictionary sd(this->subDict(typeName + "Coeffs"));

    // Get the input information and write input file
    autoPtr<OFstream> inputFile;
    inputFile.reset(new OFstream("fenton.inp"));

    inputFile().precision(14);

    // Wave height and depth information
    scalar height = readScalar(sd.lookup("height"));
    scalar depth = readScalar(sd.lookup("depth"));

    inputFile() << "finite " << height/depth << endl;

    // Specify with period scale: period/wave length
    Switch specifyPeriod(sd.lookup("specifyPeriod"));
    scalar periodScale = 0;
    if (specifyPeriod)
    {
        periodScale = readScalar(sd.lookup("period"));
        inputFile() << "period " << height/(G_*Foam::sqr(periodScale)) << endl;
    }
    else
    {
        periodScale = readScalar(sd.lookup("waveLength"));
        inputFile() << "wavelength " << height/periodScale << endl;
    }

    // Mean flow scale
    Switch specifyEuler(sd.lookup("specifyEuler"));
    scalar currentScale = 0;
    if (specifyEuler)
    {
        currentScale = readScalar(sd.lookup("eulerVelocity"));
        inputFile() << "Euler " << currentScale/Foam::sqrt(height*G_) << endl;
    }
    else
    {
        currentScale = readScalar(sd.lookup("stokesVelocity"));
        inputFile() << "Stokes " << currentScale/Foam::sqrt(height*G_) << endl;
    }

    // Get the controls
    label nModes(readLabel(sd.lookup("N")));
    label nIter(readLabel(sd.lookup("Niter")));

    N_ = nModes;

    inputFile() << nModes << " " << nIter;

    // Done writing the input file
}


void streamFunction::printCoeffs()
{
    Info << "Loading Morison force: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector streamFunction::U
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    scalar Z((x & vertDir_) - h_);
    scalar K_ = Foam::mag(k_);
    scalar Uhorz(omega_/K_ - uBar_);
    scalar Uvert(0);
    scalar arg = (k_ & x)*spaceFact_ - omega_*t;
    scalar j(0);

    forAll (B_,ii)
    {
        j = ii + 1;

        Uhorz += B_[ii]*Foam::cosh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::cos(j*arg);

        Uvert += B_[ii]*Foam::sinh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::sin(j*arg);
    }
    Uhorz *= factor;
    Uvert *= factor;

    // Note "-" because of "g" working in the opposite direction
    vector vel(Uhorz*k_/K_);

    if (includeVertical_)
    {
        vel += Uvert*vertDir_;
    }

    return dimensionedVector("null", dimVelocity, vel);
}


dimensionedVector streamFunction::acc
(
    const point& x,
    const scalar& t
) const
{
    scalar factor(rampfactor(Tsoft_));

    scalar Z((x & vertDir_) - h_);
    scalar K_ = Foam::mag(k_);
    scalar Ahorz(omega_/K_ - uBar_);
    scalar Avert(0);
    scalar arg = (k_ & x)*spaceFact_ - omega_*t;
    scalar j(0);

    forAll (B_,ii)
    {
        j = ii + 1;

        Ahorz +=  omega_*B_[ii]*Foam::cosh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::sin(j*arg);

        Avert += -omega_*B_[ii]*Foam::sinh(j*K_*(Z + h_))/Foam::cosh(j*K_*h_)
            *Foam::cos(j*arg);
    }
    Ahorz *= factor;
    Avert *= factor;

    vector acc(Ahorz*k_/K_);

    if (includeVertical_)
    {
        acc += Avert*vertDir_;
    }

    return dimensionedVector("null", dimVelocity/dimTime, acc);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
