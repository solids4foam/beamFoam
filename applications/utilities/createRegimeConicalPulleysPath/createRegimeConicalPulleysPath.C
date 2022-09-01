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
#include "HermiteSpline.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("totalPathLength");

#   include "setRootCase.H"
#   include "createTime.H"

    // Read total path length
    scalar totalPathLength
    (
        readScalar(IStringStream(args.additionalArgs()[0])())
    );

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

        if (nPulleys < 1)
        {
            FatalErrorIn(args.executable())
                << "No specified conical pulleys."
                << exit(FatalError);
        }

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
        FatalErrorIn(args.executable())
            << "No specified conical pulleys."
            << exit(FatalError);
    }

    // Create path
    DynamicField<vector> points;
    DynamicField<vector> tangents;

    //- Add first pulley points
    {
        vector axis2 =
        (
            conicalPulleys[0].axis()
           ^conicalPulleys[0].polarAxis()
        );

        
        // Calculate first point
        vector T0 =
            conicalPulleys[0].origin()
          - conicalPulleys[0].maxRadius()*axis2;
        vector tangent0 = -conicalPulleys[0].polarAxis();


        // Calculate last point
        vector T1 = T0 + 2*conicalPulleys[0].maxRadius()*axis2;
        vector tangent1 = conicalPulleys[0].polarAxis();

        
        // Add central points
        vector Tm1 =
            conicalPulleys[0].origin()
          - conicalPulleys[0].maxRadius()
           *conicalPulleys[0].polarAxis();
        
        vector tangentM1 = axis2;

        vector Tm0 = (Tm1 + T0)/2 - conicalPulleys[0].origin();
        Tm0 /= mag(Tm0) + SMALL;
        Tm0 *= conicalPulleys[0].maxRadius();
        Tm0 += conicalPulleys[0].origin();

        vector tangentM0 = (Tm1 - T0);
        tangentM0 /= mag(tangentM0) + SMALL;

        
        vector Tm2 = (T1 + Tm1)/2 - conicalPulleys[0].origin();
        Tm2 /= mag(Tm2) + SMALL;
        Tm2 *= conicalPulleys[0].maxRadius();
        Tm2 += conicalPulleys[0].origin();

        vector tangentM2 = (T1 - Tm1);
        tangentM2 /= mag(tangentM2) + SMALL;


        // Add all points
        points.append(T0);
        tangents.append(tangent0);

        points.append(Tm0);
        tangents.append(tangentM0);

        points.append(Tm1);
        tangents.append(tangentM1);

        points.append(Tm2);
        tangents.append(tangentM2);

        points.append(T1);
        tangents.append(tangent1);
    }


    //- Add second pulley points
    {
        vector axis2 =
        (
            conicalPulleys[1].axis()
           ^conicalPulleys[1].polarAxis()
        );


        // Calculate first point
        vector T2 =
            conicalPulleys[1].origin()
          + conicalPulleys[1].maxRadius()*axis2;
        vector tangent2 = conicalPulleys[1].polarAxis();


        // Calculate second point
        vector T3 = T2 - 2*conicalPulleys[1].maxRadius()*axis2;
        vector tangent3 = -conicalPulleys[1].polarAxis();

    
        // Add central points
        vector Tm1 =
            conicalPulleys[1].origin()
          + conicalPulleys[1].maxRadius()
           *conicalPulleys[1].polarAxis();

        vector tangentM1 = -axis2;

        vector Tm0 =
            (Tm1 + T2)/2
          - conicalPulleys[1].origin();
        Tm0 /= mag(Tm0) + SMALL;
        Tm0 *= conicalPulleys[1].maxRadius();
        Tm0 += conicalPulleys[1].origin();

        vector tangentM0 = (Tm1 - T2);
        tangentM0 /= mag(tangentM0) + SMALL;


        vector Tm2 =
            (T3 + Tm1)/2
          - conicalPulleys[1].origin();
        Tm2 /= mag(Tm2) + SMALL;
        Tm2 *= conicalPulleys[1].maxRadius();
        Tm2 += conicalPulleys[1].origin();

        vector tangentM2 = (T3 - Tm1);
        tangentM2 /= mag(tangentM2) + SMALL;


        // Add all points
        points.append(T2);
        tangents.append(tangent2);

        points.append(Tm0);
        tangents.append(tangentM0);

        points.append(Tm1);
        tangents.append(tangentM1);

        points.append(Tm2);
        tangents.append(tangentM2);
        
        points.append(T3);
        tangents.append(tangent3);
    }

    //Add straight part of the path at the outlet
    {
        HermiteSpline currentPath(points, tangents);
        scalar curLength = currentPath.length();
        scalar dl = totalPathLength - curLength;
        if (dl<0)
        {
            FatalErrorIn(args.executable())
                << "Total path length shuld be larger then "
                << curLength << exit(FatalError);            
        }
        
        label lastPulley = conicalPulleys.size()-1;

        vector DL =  -dl*conicalPulleys[lastPulley].polarAxis();
        DL.y() +=
            0.8
           *(
                conicalPulleys[1].maxRadius()
              - conicalPulleys[1].minRadius()
            );
        
        points.append(points[points.size()-1]+DL);
        tangents.append(tangents[tangents.size()-1]);
    }

    // Write path data to file
    {
        fileName conicalPulleysDir(runTime.caseConstant()/"conicalPulleys");
        mkDir(conicalPulleysDir);

        fileName pathFileName("conicalPulleysPath");
        OFstream pathFile(conicalPulleysDir/pathFileName);

        pathFile << points << endl << tangents << endl;
    }

    // Write path data to vtk file
    {
        HermiteSpline spline(points, tangents);

        label n = 1000;
        vectorField pts(n+1);
        scalar L = spline.length();
        scalar DL = L/n;

        Info << "L: " << L << endl;

        for (label i=0; i<=n; i++)
        {
            pts[i] = spline.position(i*DL);
        }

        fileName conicalPulleysDir(runTime.caseConstant()/"conicalPulleys");
        mkDir(conicalPulleysDir);

        fileName vtkFileName("conicalPulleysPath.vtk");
        
        OFstream vtkFile(conicalPulleysDir/vtkFileName);
        
        vtkFile << "# vtk DataFile Version 2.0" << nl
            << vtkFileName << nl
            << "ASCII" << nl
            << "DATASET UNSTRUCTURED_GRID" << nl
            // << "DATASET POLYDATA" << nl
            << "POINTS " << pts.size() << " float" << nl;
        
        forAll(pts, pointI)
        {
            vtkFile << pts[pointI].x() << ' '
                    << pts[pointI].y() << ' '
                    << pts[pointI].z() << nl;
        }

        // Write cells
        label nCells = pts.size()-1;
        vtkFile << "\nCELLS " << nCells << ' ' << nCells*3 << nl;

        for (label i=1; i<pts.size(); i++)
        {
            vtkFile << 2 << ' ' << i-1 << ' ' << i << nl;
        }

        // Write cell types
        vtkFile << "\nCELL_TYPES " << nCells << nl;

        for (label i=1; i<=nCells; i++)
        {
            vtkFile << 3 << nl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
