/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamLineContactData.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "beamModel.H"
#include "lineContact.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(beamLineContactData, 0);
    addToRunTimeSelectionTable(functionObject, beamLineContactData, dictionary);
}
}

// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

void Foam::functionObjects::beamLineContactData::makeHistoryDir()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const word startTimeName =
        time_.timeName(mesh.time().startTime().value());

    if (Pstream::parRun())
    {
        historyDir_ =
            time_.path()/".."/"postProcessing"/startTimeName/
            "beamLineContactData";
    }
    else
    {
        historyDir_ =
            time_.path()/"postProcessing"/startTimeName/"beamLineContactData";
    }

    mkDir(historyDir_);
}


bool Foam::functionObjects::beamLineContactData::writeData()
{
    const label colW = 16;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>("beamProperties");

    if (Pstream::master())
    {
        if (historyDir_.empty())
        {
            makeHistoryDir();
        }

        const fileName timeDir = historyDir_/time_.timeName();

        mkDir(timeDir);

        const auto& lineContacts = beam.contact().lineContacts();
        const label nBeams = lineContacts.size();

        forAll(lineContacts, firstBeamI)
        {
            for (label secondBeamI = 0; secondBeamI < nBeams; ++secondBeamI)
            {
                if (secondBeamI == firstBeamI)
                {
                    continue;
                }

                const word pairName =
                    "beam" + Foam::name(firstBeamI)
                  + "_to_beam" + Foam::name(secondBeamI);

                OFstream os(timeDir/(pairName + ".dat"));

                os  << setw(colW) << "segment"
                    << setw(colW) << "normalForce"
                    << setw(colW) << "forceOffset"
                    << setw(colW) << "normalGap"
                    << setw(colW) << "outerIter"
                    << endl;

                forAll(lineContacts[firstBeamI], segI)
                {
                    const lineContact& contact =
                        lineContacts[firstBeamI][segI][secondBeamI];

                    os  << setw(colW) << segI
                        << setw(colW) << mag(contact.normalContactForce())
                        << setw(colW) << contact.normalContactForceOffset()
                        << setw(colW) << contact.normalGap()
                        << setw(colW) << beam.iOuterCorr()
                        << endl;
                }
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamLineContactData::beamLineContactData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    name_(name),
    time_(runTime),
    regionName_(polyMesh::defaultRegion),
    historyDir_()
{
    Info<< "Creating " << this->name() << " function object!" << endl;

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamLineContactData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::functionObjects::beamLineContactData::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamLineContactData::end()
{
    return true;
}


bool Foam::functionObjects::beamLineContactData::write()
{
    return true;
}


// ************************************************************************* //
