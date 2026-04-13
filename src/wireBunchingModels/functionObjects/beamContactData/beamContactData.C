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

#include "beamContactData.H"
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
    defineTypeNameAndDebug(beamContactData, 0);
    addToRunTimeSelectionTable(functionObject, beamContactData, dictionary);
}
}

// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

bool Foam::functionObjects::beamContactData::writeData()
{
    const label colW = 16;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>("beamProperties");

    if (Pstream::master() && beam.contact().pointContacts().size())
    {
        const auto& pointContact = beam.contact().pointContacts()[0];

        historyFilePtr_()
        << setw(colW) << time_.time().value()
        << setw(colW) << mag(pointContact.normalContactForce())
        << setw(colW) << pointContact.normalGap()
        << setw(colW) << pointContact.normalContactForceOffset()
        << setw(colW) << pointContact.firstBeamElasticGap()
        << setw(colW) << pointContact.firstBeamSlipGap()
        << setw(colW) << pointContact.secondBeamElasticGap()
        << setw(colW) << pointContact.secondBeamSlipGap()
        << setw(colW)
        << mag
           (
               pointContact.firstBeamTangContactForce()
             - pointContact.secondBeamTangContactForce()
           )
        << setw(colW) << beam.iOuterCorr()
        << setw(colW)
        << pointContact.underRelaxationFactor()
        << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamContactData::beamContactData
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
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object!" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Create history file if not already created
    if (!historyFilePtr_.valid())
    {
        if (Pstream::master())
        {
            fileName historyDir;

            const word startTimeName =
                time_.timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                historyDir = time_.path()/".."/"postProcessing"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing"/startTimeName;
            }

            mkDir(historyDir);

            historyFilePtr_.reset
            (
                new OFstream(historyDir/"beamContactData.dat")
            );

            if (historyFilePtr_.valid())
            {
                const label colW = 16;

                historyFilePtr_()
                << setw(colW) << "Time"
                << setw(colW) << "f_con"
                << setw(colW) << "g_n"
                << setw(colW) << "f_offset"
                << setw(colW) << "fb_gE"
                << setw(colW) << "fb_gS"
                << setw(colW) << "sb_gE"
                << setw(colW) << "sb_gS"
                << setw(colW) << "net_fric"
                << setw(colW) << "nBeamCorrectors"
                << setw(colW) << "alpha"
                << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamContactData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::functionObjects::beamContactData::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamContactData::end()
{
    return true;
}


bool Foam::functionObjects::beamContactData::write()
{
    return true;
}


// ************************************************************************* //
