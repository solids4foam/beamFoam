/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

#include "beamSolverIterations.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "beamModel.H"

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

const Foam::beamModel&
Foam::functionObjects::beamSolverIterations::beam() const
{
    // Every beam model registers on Time as "beamProperties", so find the
    // one by its mesh region
    const UPtrList<const beamModel> beams
    (
        time_.csorted<beamModel>()
    );

    forAll(beams, bI)
    {
        if (beams[bI].mesh().name() == regionName_)
        {
            return beams[bI];
        }
    }

    wordList regions(beams.size());
    forAll(beams, bI)
    {
        regions[bI] = beams[bI].mesh().name();
    }

    FatalErrorInFunction
        << "No beam model on region " << regionName_ << nl
        << "Beam model regions: " << regions
        << exit(FatalError);

    return beams[0];
}


bool Foam::functionObjects::beamSolverIterations::writeData()
{
    const beamModel& b = beam();

    // Counts of a time step in which the beam was not solved are zero
    const bool solved = b.solverCountsTimeIndex() == time_.timeIndex();

    if (Pstream::master())
    {
        historyFilePtr_()
            << time_.timeName() << " "
            << (solved ? b.nEvolveCalls() : 0) << " "
            << (solved ? b.nBeamSolves() : 0) << " "
            << (solved ? b.nNewtonIterations() : 0) << " "
            << b.totalNewtonIterations()
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
            << "# Time nPimple nBeamSolves nNewton totalNewton" << endl;
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
