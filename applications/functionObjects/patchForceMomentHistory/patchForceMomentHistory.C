/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "patchForceMomentHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchForceMomentHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchForceMomentHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::patchForceMomentHistory::writeData()
{
    Info << "Writing force and moment history" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const volVectorField& W =
        mesh.lookupObject<volVectorField>("W");

    const vectorField& patchW =
        W.boundaryField()[patchIndex_];
    
    const surfaceVectorField& Q =
        mesh.lookupObject<surfaceVectorField>("Q");

    const vectorField& patchQ =
        Q.boundaryField()[patchIndex_];

    const surfaceVectorField& M =
        mesh.lookupObject<surfaceVectorField>("M");

    const vectorField& patchM =
        M.boundaryField()[patchIndex_];

    vector avgW = gAverage(patchW);
    
    vector avgQ = gSum(patchQ);
    vector avgM = gSum(patchM);
    
    Info << "Force and moment at patch " << patchName_
         << ": " << avgQ << ", " << avgM << "\n" << endl;

    if (Pstream::master())
    {
        historyFilePtr_()
            << mesh.time().value()
            << tab << avgQ.x()
            << tab << avgQ.y()
            << tab << avgQ.z()
            << tab << avgM.x()
            << tab << avgM.y()
            << tab << avgM.z()
            << tab << avgW.x()
            << tab << avgW.y()
            << tab << avgW.z() << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchForceMomentHistory::patchForceMomentHistory
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_(NULL),
    patchName_(dict.lookup("patchName")),
    patchIndex_(-1)
{
    Info << "Creating functio object " << name_ << "\n" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    polyPatchID patch(patchName_, mesh.boundaryMesh());

    if (!patch.active())
    {
        FatalErrorIn("patchForceMomentHistory::patchForceMomentHistory()")
            << "Patch name " << patchName_ << " not found."
                << abort(FatalError);
    }

    patchIndex_ = patch.index();

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                mesh.time().timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            fileName file("forceMoment_" + patchName_ + ".dat");

            // Open new file at start up
            historyFilePtr_.reset(new OFstream(historyDir/file));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time"
                    << tab << "Qx"
                    << tab << "Qy"
                    << tab << "Qz"
                    << tab << "Mx"
                    << tab << "My"
                    << tab << "Mz"
                    << tab << "Wx"
                    << tab << "Wy"
                    << tab << "Wz" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchForceMomentHistory::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::patchForceMomentHistory::execute(const bool forceWrite)
#else
bool Foam::patchForceMomentHistory::execute()
#endif
{
    return writeData();
}


bool Foam::patchForceMomentHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
