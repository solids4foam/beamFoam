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

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "faceList.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
#include "Xfer.H"
#include "IOdictionary.H"
#include "cellZone.H"

#include "conicalPulley.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"

    IOdictionary beamProperties
    (
        IOobject
        (
            "beamProperties",
            fileName(runTime.caseConstant()),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE // must be AUTO_WRITE
        )
    );

    PtrList<conicalPulley> conicalPulleys;
    
    // Read conical pulleys
    if (beamProperties.found("conicalPulleys"))
    {
        Info << "Found conical pulleys" << endl;
      
        const PtrList<entry> entries
        (
            beamProperties.lookup("conicalPulleys")
        );

        label nPulleys = entries.size();

        Info << "nPulleys: " << nPulleys << endl;
	
        conicalPulleys.setSize(nPulleys);

        forAll(entries, pulleyI)
        {
            conicalPulleys.set
            (
                pulleyI,
                new conicalPulley(entries[pulleyI].dict())
            );
        }
    }
    else
    {
    }

    fileName conicalPulleysDir(runTime.caseConstant()/"conicalPulleys");
        
    // Write OpenSCAD file
    {
        mkDir(conicalPulleysDir);

        fileName scadFileName("conicalPulleys.scad");
        
        OFstream scadFile(conicalPulleysDir/scadFileName);

        scadFile << "$fn = 100;\n" << endl;

        forAll(conicalPulleys, pulleyI)
        {
            scadFile << conicalPulleys[pulleyI];
        }
    }

    // Create stl file by calling openscad
    label rv = 0;
    {
        fileName stlFileName
        (
            runTime.caseConstant()/"conicalPulleys"/"conicalPulleys.stl"
        );

        fileName scadFileName
        (
            runTime.caseConstant()/"conicalPulleys"/"conicalPulleys.scad"
        );
        
        string command = "openscad -o " + stlFileName + " " + scadFileName;
        
        rv = system(command.c_str());
    }

    // Write OpenSCAD file for each pulley
    forAll(conicalPulleys, pulleyI)
    {
        OStringStream conicalPulleyN;
        conicalPulleyN() << "conicalPulley_" << pulleyI << ".scad";

        fileName scadFileName(conicalPulleysDir/word(conicalPulleyN.str()));
        
        OFstream scadFile(scadFileName);
        
        scadFile << "$fn = 100;\n" << endl;
        scadFile << conicalPulleys[pulleyI];

        // Create stl file
        OStringStream conicalPulleyNstl;
        conicalPulleyNstl() << "conicalPulley_" << pulleyI << ".stl";

        fileName stlFileName
        (
            runTime.caseConstant()/"conicalPulleys"
           /word(conicalPulleyNstl.str())
        );

        string command = "openscad -o " + stlFileName + " " + scadFileName;
        
        rv = system(command.c_str());
    }

    
    Info<< "End\n" << endl;

    return rv;
}


// ************************************************************************* //
