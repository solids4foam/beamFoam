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

#include "setTransverseDisplacement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "axialForceTransverseDisplacementFvPatchVectorField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setTransverseDisplacement, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setTransverseDisplacement,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setTransverseDisplacement::setBC()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    volVectorField& DW =
        const_cast<volVectorField&>
        (
            mesh.lookupObject<volVectorField>("DW")
        );

    axialForceTransverseDisplacementFvPatchVectorField& pDW =
        refCast<axialForceTransverseDisplacementFvPatchVectorField>
        (
            DW.boundaryField()[endPatchIndex_]
        );

    // scalar radius = 7.5;
    // scalar thetaEnd = 360*M_PI/180;

    scalar tEnd = time_.endTime().value();

    scalar n = tEnd/time_.deltaT().value();
    scalar dTheta = thetaEnd_/n;

    scalar z0 = -radius_;

    vector currentPosition =
        mesh.C().boundaryField()[endPatchIndex_][0];

    // Current position in polar coordinates
    scalar x = currentPosition.x();
    scalar y = currentPosition.y();
    scalar z = currentPosition.z() - z0;
    // scalar r = sqrt(sqr(y) + sqr(z));
    scalar theta = std::atan2(y, z);

    scalar newTheta = theta + dTheta;
    scalar newR = radius_;
    vector newPosition
    (
        x,
        newR*std::sin(newTheta),
        newR*std::cos(newTheta) + z0
    );

    pDW.refDisp() = newPosition - currentPosition;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setTransverseDisplacement::setTransverseDisplacement
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
    radius_(readScalar(dict.lookup("radius"))),
    thetaEnd_(readScalar(dict.lookup("angle"))*M_PI/180)
{
    Info << "Creating setTransverseDisplacement function object" << endl;
    
    if (Pstream::parRun())
    {
        FatalErrorIn("setTransverseDisplacement::setTransverseDisplacement(...)")
            << "setTransverseDisplacement objec function "
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
        FatalErrorIn("setTransverseDisplacement::setTransverseDisplacement(...)")
          << "Patch name " << startPatchName << " not found."
          << abort(FatalError);
    }
    
    startPatchIndex_ = startPatch.index();
        
    word endPatchName(dict.lookup("endPatchName"));

    polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

    if (!endPatch.active())
    {
        FatalErrorIn("setTransverseDisplacement::setTransverseDisplacement(...)")
            << "Patch name " << endPatchName << " not found."
            << abort(FatalError);
    }

    endPatchIndex_ = endPatch.index();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setTransverseDisplacement::start()
{
    return setBC();
}

bool Foam::setTransverseDisplacement::execute()
{
    return setBC();
}

bool Foam::setTransverseDisplacement::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
