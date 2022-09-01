/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Class
    plasticityStressResultantReturn

\*---------------------------------------------------------------------------*/

#include "plasticityStressResultantReturn.H"
#include "volFields.H"
#include "fvc.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(plasticityStressResultantReturn, 0);
defineRunTimeSelectionTable(plasticityStressResultantReturn, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasticityStressResultantReturn::plasticityStressResultantReturn
(
    const word& name,
    beamModel& beamModel
)
:
    name_(name),
    plasticityProperties_(beamModel.beamProperties().subDict(name + "Coeffs")),
    curTimeIndex_(-1),
    nElPredCorr_(plasticityProperties_.lookupOrDefault<int>("nElPredCorr", 10)),
    elPredTol_(plasticityProperties_.lookupOrDefault<scalar>("elPredTol", 1e-4))
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
