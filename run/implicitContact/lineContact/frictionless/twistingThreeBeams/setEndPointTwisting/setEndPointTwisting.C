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

#include "setEndPointTwisting.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "axialForceTransverseDisplacementNRFvPatchVectorField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setEndPointTwisting, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setEndPointTwisting,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setEndPointTwisting::setBC()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    scalar tEnd = time_.endTime().value();

    vectorField curPosition =
        mesh.C().boundaryField()[endPatchIndex_];

    vector i(1, 0, 0);
    scalarField curRadius = mag((I-i*i) & curPosition);
    Info << regionName_ << ", curRadius: " << curRadius << endl;
        
    {
        volVectorField& W =
            const_cast<volVectorField&>
            (
                mesh.lookupObject<volVectorField>("W")
            );

        scalar Theta =
            thetaEnd_*((time_.value() + time_.deltaT().value())/tEnd);

        if (mesh.moving()) // Updated Lagrangian formulation
        {
            scalar oldTheta =
                thetaEnd_*(time_.value()/tEnd);

            Theta -= oldTheta;
        }
        
        Info << "Theta: " << Theta << endl;
        
        // Rotation tensor
        tensor A
        (
            1, 0, 0,
            0, ::cos(Theta), -::sin(Theta),
            0, ::sin(Theta),  ::cos(Theta)
        );

        // New end point position
        vectorField newPosition = (A & curPosition);
        newPosition.replace(0, 5.049647);

        if
        (
            isA<axialForceTransverseDisplacementNRFvPatchVectorField>
            (
                W.boundaryField()[endPatchIndex_]
            )
        )
        {
            axialForceTransverseDisplacementNRFvPatchVectorField& pW =
                refCast<axialForceTransverseDisplacementNRFvPatchVectorField>
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
    }
    
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setEndPointTwisting::setEndPointTwisting
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
    Info << "Creating setEndPointTwisting function object" << endl;
    
    if (Pstream::parRun())
    {
        FatalErrorIn("setEndPointTwisting::setEndPointTwisting(...)")
            << "setEndPointTwisting objec function "
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
        FatalErrorIn("setEndPointTwisting::setEndPointTwisting(...)")
          << "Patch name " << startPatchName << " not found."
          << abort(FatalError);
    }

    startPatchIndex_ = startPatch.index();

    word endPatchName(dict.lookup("endPatchName"));

    polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

    if (!endPatch.active())
    {
        FatalErrorIn("setEndPointTwisting::setEndPointTwisting(...)")
            << "Patch name " << endPatchName << " not found."
            << abort(FatalError);
    }

    endPatchIndex_ = endPatch.index();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setEndPointTwisting::start()
{
    return setBC();
}

#if FOAMEXTEND > 40
bool Foam::setEndPointTwisting::execute(const bool forceWrite)
#else
bool Foam::setEndPointTwisting::execute()
#endif
{
    return setBC();
}

bool Foam::setEndPointTwisting::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
