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
    Path calculation is performed using procedure described in:
    http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm

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
#include "DynamicField.H"
#include "cylindricalCS.H"

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

    // Create path
    DynamicField<vector> points;
    DynamicField<vector> tangents;

    //- Add first point 
    vector axis2 =
    (
        conicalPulleys[0].polarAxis()
       ^conicalPulleys[0].axis()
    );
    points.append
    (
        conicalPulleys[0].origin()
      + conicalPulleys[0].maxRadius()*axis2
    );
    tangents.append(conicalPulleys[0].polarAxis());

    //- Add remining points
    for (label i=1; i < conicalPulleys.size(); i++)
    {
        cylindricalCS cs
        (
            "cs",
            conicalPulleys[i-1].origin(),
            conicalPulleys[i-1].axis(),
            conicalPulleys[i-1].polarAxis()
        );

        vector C0 = cs.globalToLocal(conicalPulleys[i-1].origin());
        vector C1 = cs.globalToLocal(conicalPulleys[i].origin());
        
        scalar a = C0.x();
        scalar b = C0.y();

        scalar c = C1.x();
        scalar d = C2.y();

        scalar r0 = conicalPulleys[i-1].maxRadius();
        scalar r1 = conicalPulleys[i].maxRadius();
        
        scalar xP = (c*r0 + a*r1)/(r0+r1);
        scalar yP = (d*r0 + b*r1)/(r0+r1);

        scalar xt1
    }
    
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
