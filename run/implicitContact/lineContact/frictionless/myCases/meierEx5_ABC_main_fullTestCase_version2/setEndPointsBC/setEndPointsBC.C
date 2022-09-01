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

#include "setEndPointsBC.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "axialForceTransverseDisplacementFvPatchVectorField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setEndPointsBC, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setEndPointsBC,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setEndPointsBC::setBC()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    scalar tEnd = time_.endTime().value();
    
    label nBeams = mesh.cellZones().size();
    labelList nZ;
    nZ.setSize(nBeams);
    
    for (label bI = 0; bI < nBeams; bI++)
    {
	nZ[bI] = mesh.cellZones()[bI].size();
    }
    
    Info << "nBeams " << nBeams << " nZ " << nZ[0]+0.5*(nZ[1]-1) << endl;

    //   vectorField curPosition =
        //   mesh.C().boundaryField()[endPatchIndex_];
    {
        Info << "Set BC for W and Theta in total Lagrangian solver" << endl;
	
	volVectorField& W =
	const_cast<volVectorField&>
	(
	    mesh.lookupObject<volVectorField>("W")
	);
		
	if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                W.boundaryField()[endPatchIndex_]
            )
        )
        {	    
	    FatalErrorIn("setEndPointsBC::setBC(...) ")
		<< "setEndPointsBC function object "
		<< " is not implemented for "
		<< " axialForceFVpatchVetor"
		<< abort(FatalError);
        }
        else
        {
	    volVectorField& Theta =
	    const_cast<volVectorField&>
	    (
		mesh.lookupObject<volVectorField>("Theta")
	    );
	    
	    
	    if (mesh.moving()) // Updated Lagrangian formulation
	    {

		FatalErrorIn("setEndPointsBC::setBC(...) ")
		<< "setEndPointsBC function object "
		<< " is not implemented for "
		<< " Updated Lagrangian formulation"
		<< abort(FatalError);
   
	    }
	
	    if (time_.value() < 1)
	    {
		Info << "time " << time_.value() << endl;
		
		scalar curYdispMag =
		    yDispEnd_*(time_.value() + time_.deltaT().value());
		
		Info << "current downward disp " << curYdispMag << endl;
		    
		vector curTheta(0, 0, 0);

		Theta.internalField()[nZ[0]+0.5*(nZ[1] - 1)] == curTheta;
     	
		vector curW(0, curYdispMag, 0);
		
		W.boundaryField()[startPatchIndex_] == curW;
		W.boundaryField()[endPatchIndex_] == curW; 

		
	    }
	    else
	    {
		Info << "Inside time > 1 sec loop" << endl;
		
		scalar thetaMag =
		    (thetaEnd_/(tEnd - 1))*(time_.value() - 1);
		    	    
		// Setting displacement at end patches
		vector curWleft
		(
		-::sin(thetaMag), yDispEnd_, -(1 - ::cos(thetaMag))
		);
		
		vector curWright
		(
		::sin(thetaMag), yDispEnd_, (1 - ::cos(thetaMag))
		);

		W.boundaryField()[startPatchIndex_] == curWleft;
		W.boundaryField()[endPatchIndex_] == curWright; 
		
		// Setting the rotation
		vector curTheta(0, thetaMag, 0);
	
		Theta.boundaryField()[startPatchIndex_] == -curTheta;
	    }
		    		
        }
	    
	//   scalar curYdispMag = 
	       //   yDispEnd_*((time_.value() + time_.deltaT().value())/tEnd);

	//   vector curTheta(-::asin(curYdispMag), 0, 0);
	//   vector curTheta(0, 0, 0);
	
	//   Info << "Theta middle before " << Theta.internalField()[150]
	     //   << endl;
	     
	//   if (curYdispMag < -0.1)
	//   {
	    //   vector curTheta(-::asin(curYdispMag), 0, 0);
	    //   Theta.boundaryField()[startPatchIndex_] == curTheta;
	    //   Theta.boundaryField()[endPatchIndex_] == curTheta;
	//   }
	//   else
	//   {
	    //   vector curTheta(0, 0, 0);
	 //   Theta.internalField()[150] == curTheta;
	//   }
	
	
	//Theta.boundaryField()[endPatchIndex_] == -curTheta;
    }
    
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setEndPointsBC::setEndPointsBC
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
    //radius_(-1),
    // radius_(readScalar(dict.lookup("radius"))),
    yDispEnd_(readScalar(dict.lookup("yDisp"))),
    // thetaEnd_(readScalar(dict.lookup("angle"))*M_PI/180)
{
    Info << "Creating setEndPointsBC function object" << endl;

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
        FatalErrorIn("setEndPointsBC::setEndPointsBC(...)")
          << "Patch name " << startPatchName << " not found."
          << abort(FatalError);
    }

    startPatchIndex_ = startPatch.index();

    word endPatchName(dict.lookup("endPatchName"));

    polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

    if (!endPatch.active())
    {
        FatalErrorIn("setEndPointsBC::setEndPointsBC(...)")
            << "Patch name " << endPatchName << " not found."
            << abort(FatalError);
    }

    endPatchIndex_ = endPatch.index();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setEndPointsBC::start()
{
    return setBC();
}

#if FOAMEXTEND > 40
bool Foam::setEndPointsBC::execute(const bool forceWrite)
#else
bool Foam::setEndPointsBC::execute()
#endif
{
    return setBC();
}

bool Foam::setEndPointsBC::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
