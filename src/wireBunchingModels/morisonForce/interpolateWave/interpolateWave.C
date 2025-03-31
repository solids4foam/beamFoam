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

#include "interpolateWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(interpolateWave, 0);
addToRunTimeSelectionTable(morisonForce, interpolateWave, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


interpolateWave::interpolateWave
(
    const Time& runTime,
    const fvMesh& mesh,
    const dictionary& sd
)
:
    morisonForce(runTime, mesh, sd),

    z_(this->subDict(typeName + "Coeffs").lookup("z")),

    t_(this->subDict(typeName + "Coeffs").lookup("time")),

    invc_(this->subDict(typeName + "Coeffs").lookup("waveNumber")),

    velMean_(this->subDict(typeName + "Coeffs").lookupOrDefault<vector>("velMean", vector::zero))
{
    // Define the size of the velocity and acceleration fields
    velocity_.setSize(z_.size());
    acceleration_.setSize(z_.size());

    // Loop over all coordinates and read data
    forAll (z_, zI)
    {
        // Read the velocity field
        velocity_[zI].setSize(t_.size());

        word name("velocity" + std::to_string(zI));

        vectorField tmpV(name, this->subDict(typeName + "Coeffs"), t_.size());
        velocity_[zI] = tmpV;

        // Read the acceleration field
        acceleration_[zI].setSize(t_.size());

        word nameA("acceleration" + std::to_string(zI));

        vectorField tmpA(nameA, this->subDict(typeName + "Coeffs"), t_.size());
        acceleration_[zI] = tmpA;
    }


    // Make the invC vector
    invc_ /= readScalar(this->subDict(typeName + "Coeffs").lookup("omega"));
}


void interpolateWave::printCoeffs()
{
    Info << "Loading Morison force: " << typeName << endl;
}


vector interpolateWave::interpolateField
(
    const List<vectorField>& field,
    const scalar& z,
    const scalar t
) const
{
    // Identify the z-interval
    label zInt(0);

    for (int ii=0; ii < z_.size() - 1; ii++)
    {
        if (z_[ii] <= z && z <= z_[ii+1])
        {
            zInt = ii;
        }
    }

    // Identify the t-time under the assumption that the time axis is equidistant
    scalar dt(t_[1] - t_[0]);
    long Nt(std::floor((t - t_[0])/dt));
    scalar w((t_[Nt + 1] - t)/dt);

    // Calculate the field at the two levels of the z-interval
    vector field0(w*field[zInt][Nt] + (1. - w)*field[zInt][Nt + 1]);
    vector field1(w*field[zInt + 1][Nt] + (1. - w)*field[zInt + 1][Nt + 1]);

    // Calculate the vertical interpolation weight
    w = (z_[zInt + 1] - z)/(z_[zInt + 1] - z_[zInt]);

    return w*field0 + (1 - w)*field1;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

dimensionedVector interpolateWave::U
(
    const point& x,
    const scalar& t
)const
{
    scalar z(x & vertDir_);

    scalar tEff(t - (invc_ & x));

    return dimensionedVector
        (
            "null",
            dimVelocity,
            this->interpolateField(velocity_, z, tEff) + velMean_
        );
}


dimensionedVector interpolateWave::acc
(
    const point& x,
    const scalar& t
) const
{
    scalar z(x & vertDir_);

    scalar tEff(t - (invc_ & x));

    return dimensionedVector
        (
            "null",
            dimVelocity,
            this->interpolateField(acceleration_, z, tEff)
        );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
