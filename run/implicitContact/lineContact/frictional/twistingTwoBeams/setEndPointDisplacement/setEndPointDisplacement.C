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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "setEndPointDisplacement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "axialForceTransverseDisplacementFvPatchVectorField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setEndPointDisplacement, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setEndPointDisplacement,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setEndPointDisplacement::setBC()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    scalar tEnd = time_.endTime().value() - time_.deltaT().value();

    vectorField curPosition =
        mesh.C().boundaryField()[endPatchIndex_];

    vector i(1, 0, 0);
    scalarField curRadius = mag((I-i*i) & curPosition);
    Info << regionName_ << ", curRadius: " << curRadius << endl;
        
    volVectorField& W =
        const_cast<volVectorField&>
        (
	    mesh.lookupObject<volVectorField>("W")
	);

    // Info << "Time index: " <<  time_.timeIndex() << endl;
    
    if (time_.timeIndex() == 0)
    {
        Info << "Set axial load in first time step" << endl;
        W.boundaryField()[endPatchIndex_] == vector(0.01, 0, 0);
    }
    else
    {
        scalar Theta = thetaEnd_*time_.value()/tEnd;
            // thetaEnd_*((time_.value() + time_.deltaT().value())/tEnd);

        // Info << "Theta0: " << Theta << endl;
	
        if (mesh.moving()) // Updated Lagrangian formulation
        {
            Info << "Set BC for W in updated Lagrangian solver" << endl;
            scalar oldTheta =
                thetaEnd_*((time_.value()-time_.deltaT().value())/tEnd);

            // Info << "oldTheta: " << oldTheta << endl;
	
            Theta -= oldTheta;
        }
	else
	{
            Info << "Set BC for W in total Lagrangian solver" << endl;
	}
        // Info << "Theta: " << Theta << endl;
        
        // Rotation tensor
        tensor A
        (
            1, 0, 0,
            0, ::cos(Theta), -::sin(Theta),
            0, ::sin(Theta),  ::cos(Theta)
        );

        // New end point position
        vectorField newPosition = (A & curPosition);
        newPosition.replace(0, 5.01);

        if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                W.boundaryField()[endPatchIndex_]
            )
        )
        {
            axialForceTransverseDisplacementFvPatchVectorField& pW =
                refCast<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    W.boundaryField()[endPatchIndex_]
                );

            pW.refDisp() = newPosition - curPosition;
        }
        else
        {
            W.boundaryField()[endPatchIndex_] ==
                (newPosition - curPosition);
        }

        // Info << (newPosition - curPosition) << endl;
    }

    // Write torque history
    {
        const surfaceVectorField& Q =
            mesh.lookupObject<surfaceVectorField>("Q");

	vector M = gSum(curPosition ^ Q.boundaryField()[endPatchIndex_]);

        if (Pstream::master())
        {
	    historyFilePtr_()
	        << time_.time().value() << " "
		<< M.x() << endl;
        }
    }
    
    
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setEndPointDisplacement::setEndPointDisplacement
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
    startPatchIndex_(-1),
    endPatchIndex_(-1),
    radius_(-1),
    // radius_(readScalar(dict.lookup("radius"))),
    thetaEnd_(readScalar(dict.lookup("angle"))*M_PI/180)
{
    Info << "Creating setEndPointDisplacement function object" << endl;
    
    if (Pstream::parRun())
    {
        FatalErrorIn("setEndPointDisplacement::setEndPointDisplacement(...)")
            << "setEndPointDisplacement objec function "
                << "is not implemented for parallel run"
                << abort(FatalError);
    }

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    word startPatchName(dict.lookup("startPatchName"));

    polyPatchID startPatch(startPatchName, mesh.boundaryMesh());
    
    if (!startPatch.active())
    {
        FatalErrorIn("setEndPointDisplacement::setEndPointDisplacement(...)")
          << "Patch name " << startPatchName << " not found."
          << abort(FatalError);
    }

    startPatchIndex_ = startPatch.index();

    word endPatchName(dict.lookup("endPatchName"));

    polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

    if (!endPatch.active())
    {
        FatalErrorIn("setEndPointDisplacement::setEndPointDisplacement(...)")
            << "Patch name " << endPatchName << " not found."
            << abort(FatalError);
    }

    endPatchIndex_ = endPatch.index();
    
    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());


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

            // Open new file at start up
            historyFilePtr_.reset
                (
                    new OFstream
                    (
                        historyDir/"endForceMoment_" + endPatchName + ".dat"
                    )
                );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time Moment" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setEndPointDisplacement::start()
{
    return setBC();
}

bool Foam::setEndPointDisplacement::execute()
{
    return setBC();
}

bool Foam::setEndPointDisplacement::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
