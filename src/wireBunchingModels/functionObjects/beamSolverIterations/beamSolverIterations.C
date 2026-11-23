/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamSolverIterations.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "beamModel.H"
#include "lookupBeamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(beamSolverIterations, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        beamSolverIterations,
        dictionary
    );
}
}

// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

bool Foam::functionObjects::beamSolverIterations::writeData()
{
    const beamSolverCounters& counts =
        lookupBeamModel(time_.lookupObject<fvMesh>(regionName_))
       .solverCounters();

    const label timeIndex = time_.timeIndex();

    if (Pstream::master())
    {
        historyFilePtr_()
            << time_.timeName() << " "
            << counts.nEvolveCalls(timeIndex) << " "
            << counts.nNewtonIterations(timeIndex) << " "
            << counts.totalNewtonIterations()
            << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamSolverIterations::beamSolverIterations
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    time_(runTime),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_()
{
    read(dict);

    if (Pstream::master())
    {
        fileName historyDir;

        const word startTimeName =
            time_.timeName(time_.startTime().value());

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            historyDir = time_.path()/".."/"postProcessing"/startTimeName;
        }
        else
        {
            historyDir = time_.path()/"postProcessing"/startTimeName;
        }

        mkDir(historyDir);

        historyFilePtr_.reset(new OFstream(historyDir/(name + ".dat")));

        historyFilePtr_()
            << "# Time nEvolve nNewton totalNewton" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamSolverIterations::read(const dictionary& dict)
{
    dict.readIfPresent("region", regionName_);

    return true;
}


bool Foam::functionObjects::beamSolverIterations::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamSolverIterations::write()
{
    return true;
}


// ************************************************************************* //
