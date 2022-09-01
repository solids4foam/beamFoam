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

#include "setEndPointsDisplacement.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "axialForceTransverseDisplacementFvPatchVectorField.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setEndPointsDisplacement, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setEndPointsDisplacement,
        dictionary
    );
    

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setEndPointsDisplacement::setBC()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);
	
    //   IOdictionary beamProperties
    //   (
        //   IOobject
        //   (
            //   "beamProperties",
            //   time_.caseConstant(),
            //   IOobject::MUST_READ,
            //   IOobject::NO_WRITE
        //   )
    //   );
    
	
    label nBeams = mesh.cellZones().size();
    
    //   scalar nSeg = mesh.cellZones()[1].size();
    
    
    labelList nZ;
    nZ.setSize(nBeams);
    
    for (label bI = 0; bI < nBeams; bI++)
    {
	nZ[bI] = mesh.cellZones()[bI].size();
    }
    
     Info << "nBeams " << nBeams << " nZ " << nZ[0]+0.5*(nZ[1]-1) << endl;

    scalar tEnd = time_.endTime().value();

    //   vectorField curPosition =
        //   mesh.C().boundaryField()[endPatchIndex_];
    {
        Info << "Set BC for W and Theta in total Lagrangian solver" << endl;
        
        volVectorField& W =
            const_cast<volVectorField&>
            (
                mesh.lookupObject<volVectorField>("W")
            );
	    
	volVectorField& Theta =
	    const_cast<volVectorField&>
	    (
		mesh.lookupObject<volVectorField>("Theta")
	    );

       // scalar thetaMag =
        //    thetaEnd_*((time_.value() + time_.deltaT().value())/tEnd);
	    
	   scalar curYdispMag = 
	       yDispEnd_*((time_.value() + time_.deltaT().value())/tEnd);

        if (mesh.moving()) // Updated Lagrangian formulation
        {

	FatalErrorIn("setEndPointsDisplacement::setDisplacement(...) ")
          << "setEndPointsDisplacement function object "
              << " is not implemented for "
	      << " Updated Lagrangian formulation"
              << abort(FatalError);
   
        }
        
        // Pout << "Theta: " << Theta << endl;
	
	//   scalarField thetax = Theta.component(vector::X);
	//   scalarField thetay = Theta.component(vector::Y);
	//   scalarField thetaz = Theta.component(vector::Z);
	
	//   vector curTheta(-::asin(curYdispMag), 0, 0);
	vector curTheta(0, 0, 0);
	
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
	Theta.internalField()[nZ[0]+0.5*(nZ[1] - 1)] == curTheta;
	//   Theta.internalField()[nZ[0]-1+0.5*(nZ[1] - 1)] == curTheta;
	//   }
	
	
	//Theta.boundaryField()[endPatchIndex_] == -curTheta;
	

	     	
	vector curW(0, curYdispMag, 0);
     
        if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                W.boundaryField()[endPatchIndex_]
            )
        )
        {	    
	    FatalErrorIn("setEndPointsDisplacement::setBC(...) ")
		<< "setEndPointsDisplacement function object "
		<< " is not implemented for "
		<< " axialForceFVpatchVetor"
		<< abort(FatalError);
        }
        else
        {	
	    W.boundaryField()[startPatchIndex_] == curW;
	    W.boundaryField()[endPatchIndex_] == curW; 
        }
    }
    
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setEndPointsDisplacement::setEndPointsDisplacement
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
    yDispEnd_(readScalar(dict.lookup("yDisp")))
    // thetaEnd_(readScalar(dict.lookup("angle"))*M_PI/180)
{
    Info << "Creating setEndPointsDisplacement function object" << endl;

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
        FatalErrorIn("setEndPointsDisplacement::setEndPointsDisplacement(...)")
          << "Patch name " << startPatchName << " not found."
          << abort(FatalError);
    }

    startPatchIndex_ = startPatch.index();

    word endPatchName(dict.lookup("endPatchName"));

    polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

    if (!endPatch.active())
    {
        FatalErrorIn("setEndPointsDisplacement::setEndPointsDisplacement(...)")
            << "Patch name " << endPatchName << " not found."
            << abort(FatalError);
    }

    endPatchIndex_ = endPatch.index();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setEndPointsDisplacement::start()
{
    return setBC();
}

#if FOAMEXTEND > 40
bool Foam::setEndPointsDisplacement::execute(const bool forceWrite)
#else
bool Foam::setEndPointsDisplacement::execute()
#endif
{
    return setBC();
}

bool Foam::setEndPointsDisplacement::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
