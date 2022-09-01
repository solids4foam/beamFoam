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

#include "toroidalPulley.H"
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

    PtrList<toroidalPulley> toroidalPulleys;

    // Read toroidal pulleys
    if (beamProperties.found("toroidalPulleys"))
    {
        Info << "Found toroidal pulleys" << endl;
      
        const PtrList<entry> entries
        (
            beamProperties.lookup("toroidalPulleys")
        );

        label nPulleys = entries.size();

        Info << "nPulleys: " << nPulleys << endl;
	
        toroidalPulleys.setSize(nPulleys);

        forAll(entries, pulleyI)
        {
            toroidalPulleys.set
            (
                pulleyI,
                new toroidalPulley(entries[pulleyI].dict())
            );
        }
    }
    else
    {
    }

    fileName toroidalPulleysDir(runTime.caseConstant()/"toroidalPulleys");
    mkDir(toroidalPulleysDir);
        
    // Write OpenSCAD file
    {
        fileName scadFileName("toroidalPulleys.scad");
        
        OFstream scadFile(toroidalPulleysDir/scadFileName);

        scadFile << "$fn = 150;\n" << endl;

        forAll(toroidalPulleys, pulleyI)
        {
            scadFile << toroidalPulleys[pulleyI];
        }
    }

    // Create stl file by calling openscad
    label rv = 0;
    {
        fileName stlFileName
        (
            runTime.caseConstant()/"toroidalPulleys"/"toroidalPulleys.stl"
        );
        
        fileName scadFileName
        (
            runTime.caseConstant()/"toroidalPulleys"/"toroidalPulleys.scad"
        );
        
        string command = "openscad -o " + stlFileName + " " + scadFileName;
        
        rv = system(command.c_str());
    }
    
    // Write OpenSCAD file for each pulley
    forAll(toroidalPulleys, pulleyI)
    {
        OStringStream toroidalPulleyN;
        toroidalPulleyN() << "toroidalPulley_" << pulleyI << ".scad";

        fileName scadFileName(toroidalPulleysDir/word(toroidalPulleyN.str()));
        
        OFstream scadFile(scadFileName);
        
        scadFile << "$fn = 150;\n" << endl;
        scadFile << toroidalPulleys[pulleyI];

        // Create stl file
        OStringStream toroidalPulleyNstl;
        toroidalPulleyNstl() << "toroidalPulley_" << pulleyI << ".stl";

        fileName stlFileName
        (
            runTime.caseConstant()/"toroidalPulleys"
           /word(toroidalPulleyNstl.str())
        );

        string command = "openscad -o " + stlFileName + " " + scadFileName;

        // Info << command << endl;
        
        rv = system(command.c_str());
    }
    
    Info<< "End\n" << endl;

    return rv;
}


// ************************************************************************* //
