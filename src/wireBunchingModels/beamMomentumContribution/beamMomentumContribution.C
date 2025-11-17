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

#include "beamMomentumContribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(beamMomentumContribution, 0);
defineRunTimeSelectionTable(beamMomentumContribution, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


autoPtr<beamMomentumContribution> beamMomentumContribution::New
(
    // const Time& runTime
    const word& beamMomentumContributionTypeName,
    const dictionary& dict
)
{
    // Info<< "Beam momentum contribution type: "
    //     << beamMomentumContributionTypeName << endl;

    // Enclose the creation of the dict to ensure it is
    // deleted before the beamMomentumContribution is created otherwise the dictionary
    // is entered in the database twice
    // {
    //     IOdictionary dict
    //     (
    //         IOobject
    //         (
    //             "beamMomentumContributionProperties",
    //             runTime.constant(),
    //             runTime,
    //             IOobject::MUST_READ,
    //             IOobject::NO_WRITE
    //         )
    //     );

    //     dict.lookup("type")
    //         >> beamMomentumContributionTypeName;
    // }

    auto* cstrIter = dictionaryConstructorTable
    (
        beamMomentumContributionTypeName
    );

    if (!cstrIter)
    {
        FatalIOErrorIn
        (
            "beamMomentumContribution::New(\n"
            "    const word& name,\n"
            "    const dictionary& dict\n"
            ")",
            dict
        )   << "Unknown beamMomentumContribution type "
            << beamMomentumContributionTypeName
            << endl << endl
            << "Valid beamMomentumContribution types are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    // auto* ctorPtr = dictionaryConstructorTable(beamMomentumContributionTypeName);

    // if (!ctorPtr)
    // {
    //     FatalErrorInFunction
    //         << "Cannot find type = " << beamMomentumContributionTypeName
    //         << exit(FatalIOError);
    // }

    return autoPtr<beamMomentumContribution>
    (
        cstrIter
            (
                 beamMomentumContributionTypeName,
                 dict
            )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


beamMomentumContribution::beamMomentumContribution
(
    const word& name,
    const dictionary& dict
)
:
//     IOdictionary dict
//     (
//         IOobject
//         (
//             "beamMomentumContributionProperties",
//             runTime.constant(),
//             runTime,
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//      ),
    name_(name)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //


beamMomentumContribution::~beamMomentumContribution()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
