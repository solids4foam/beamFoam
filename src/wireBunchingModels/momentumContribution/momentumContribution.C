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

#include "momentumContribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(momentumContribution, 0);
defineRunTimeSelectionTable(momentumContribution, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<momentumContribution> momentumContribution::New
(
    const Time& runTime
)
{
    word momentumContributionTypeName;

    // Enclose the creation of the dict to ensure it is
    // deleted before the momentumContribution is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary dict
        (
            IOobject
            (
                "momentumContributionProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        dict.lookup("type")
            >> momentumContributionTypeName;
    }

    auto* ctorPtr = dictionaryConstructorTable(momentumContributionTypeName);

    if (!ctorPtr)
    {
        FatalErrorInFunction
            << "Cannot find type = " << momentumContributionTypeName
            << exit(FatalIOError);
    }

    return autoPtr<momentumContribution>(ctorPtr(runTime));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


momentumContribution::momentumContribution
(
    const Time& runTime
)
:
    IOdictionary
    (
        IOobject
        (
            "momentumContributionProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //


momentumContribution::~momentumContribution()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
