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

#include "beamPointContactData.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(beamPointContactData, 0);
    addToRunTimeSelectionTable(functionObject, beamPointContactData, dictionary);
}
}

// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

void Foam::functionObjects::beamPointContactData::makeHistoryDir()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const word startTimeName =
        time_.timeName(mesh.time().startTime().value());

    if (Pstream::parRun())
    {
        historyDir_ = time_.path()/".."/"postProcessing"/startTimeName;
    }
    else
    {
        historyDir_ = time_.path()/"postProcessing"/startTimeName;
    }

    mkDir(historyDir_);
}


void Foam::functionObjects::beamPointContactData::makeHistoryFile
(
    const label pointContactI
)
{
    if (!Pstream::master())
    {
        return;
    }

    if (historyDir_.empty())
    {
        makeHistoryDir();
    }

    if (pointContactI >= historyFilePtrs_.size())
    {
        historyFilePtrs_.setSize(pointContactI + 1);
    }

    if (historyFilePtrs_.set(pointContactI))
    {
        return;
    }

    historyFilePtrs_.set
    (
        pointContactI,
        new OFstream
        (
            historyDir_/
            (
                "beamPointContactData_"
              + Foam::name(pointContactI)
              + ".dat"
            )
        )
    );

    if (historyFilePtrs_.set(pointContactI))
    {
        const label colW = 16;

        historyFilePtrs_[pointContactI]
            << setw(colW) << "Time"
            << setw(colW) << "normalForce"
            << setw(colW) << "forceOffset"
            << setw(colW) << "normalGap"
            << setw(colW) << "outerIter"
            << endl;
    }
}


bool Foam::functionObjects::beamPointContactData::writeData()
{
    const label colW = 16;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>("beamProperties");

    if (Pstream::master())
    {
        const auto& pointContacts = beam.contact().pointContacts();

        forAll(pointContacts, pointContactI)
        {
            makeHistoryFile(pointContactI);

            const auto& pointContact = pointContacts[pointContactI];

            historyFilePtrs_[pointContactI]
                << setw(colW) << time_.time().value()
                << setw(colW) << mag(pointContact.normalContactForce())
                << setw(colW) << pointContact.normalContactForceOffset()
                << setw(colW) << pointContact.normalGap()
                << setw(colW) << beam.iOuterCorr()
                << endl;
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamPointContactData::beamPointContactData
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
    historyDir_(),
    historyFilePtrs_()
{
    Info<< "Creating " << this->name() << " function object!" << endl;

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamPointContactData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::functionObjects::beamPointContactData::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamPointContactData::end()
{
    return true;
}


bool Foam::functionObjects::beamPointContactData::write()
{
    return true;
}


// ************************************************************************* //
