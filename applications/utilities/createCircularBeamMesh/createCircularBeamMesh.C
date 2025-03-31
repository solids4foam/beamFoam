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
Custom beam meshing utility 

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "faceList.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
//#include "Xfer.H"
#include "IOdictionary.H"
#include "cellZone.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    word regionName;
    fileName polyMeshDir;

    if (args.found("region"))
    {
        // constant/<region>/polyMesh/blockMeshDict
        regionName  = args.get("region");
        polyMeshDir = regionName/polyMesh::meshSubDir;

        Info<< nl << "Generating mesh for region " << regionName << endl;
    }
    else
    {
        // constant/polyMesh/blockMeshDict
        regionName  = polyMesh::defaultRegion;
        polyMeshDir = polyMesh::meshSubDir;
    }

    IOdictionary beamProperties
    (
        IOobject
        (
            "beamProperties",
            bool(regionName == polyMesh::defaultRegion)
          ? fileName(runTime.caseConstant())
          : fileName(runTime.caseConstant()/regionName),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE // must be AUTO_WRITE
        )
    );

    const scalar R = dimensionedScalar("R", beamProperties).value();
    const scalar L = dimensionedScalar("L", beamProperties).value();
    // scalar R = 0.1;
    // scalar L = 20;

    const label nZ(readInt(beamProperties.lookup("nSegments")));
    // label nZ = 50;
    const label nTheta = 32;

    const label nPoints = nTheta*(nZ + 1);
    const label nFaces = (nZ+1) + nZ*nTheta;
    const label nInternalFaces = nZ - 1;


    // Points
    // Xfer<vectorField> points(vectorField(nPoints, vector::zero));

    vectorField points(nPoints, vector::zero);


    scalar dTheta = -2*M_PI/nTheta;

    vectorField circumferencePoints(nTheta, vector::zero);
    forAll(circumferencePoints, pointI)
    {
        circumferencePoints[pointI].x() = R*::cos(pointI*dTheta);
        circumferencePoints[pointI].y() = R*::sin(pointI*dTheta);
        circumferencePoints[pointI].z() = 0;
    }

    label pI = 0;
    scalar dZ = L/nZ;
    for (label i=0; i<=nZ; i++)
    {
        for (label j=0; j<nTheta; j++)
        {
            vector curPoint = circumferencePoints[j];
            curPoint.z() = i*dZ;

            points[pI++] = curPoint;
        }
    }


    // Faces
    faceList faces(nFaces, face());

    //-Internal faces

    label fI = 0;
    for (label i=1; i<=nInternalFaces; i++)
    {
        labelList curPointLabels(nTheta, -1);
        for (label j=0; j<nTheta; j++)
        {
            curPointLabels[j] = j + i*nTheta;
        }

        face curFace(curPointLabels);
        curFace = curFace.reverseFace();
        //faces()[fI++] = curFace;
	    faces[fI++] = curFace;
        // faces()[fI++] = face(curPointLabels);
    }

    //- Inlet face

    labelList curPointLabels(nTheta, -1);
    for (label j=0; j<nTheta; j++)
    {
        curPointLabels[j] = j;
    }

    face curFace(curPointLabels);
    // curFace = curFace.reverseFace();

    faces[fI++] = curFace;

    //- Outlet face

    curPointLabels = labelList(nTheta, -1);
    for (label j=0; j<nTheta; j++)
    {
        curPointLabels[j] = j + nZ*nTheta;
    }

    curFace = face(curPointLabels);
    curFace = curFace.reverseFace();
    faces[fI++] = curFace;

    // faces()[fI++] = face(curPointLabels);

    //- Pipe faces

    for (label i=0; i<nZ; i++)
    {
        for (label j=1; j<nTheta; j++)
        {
            labelList curPointLabels(4, -1);

            curPointLabels[0] = (j-1) + i*nTheta;
            curPointLabels[1] = j + i*nTheta;
            curPointLabels[2] = j + (i+1)*nTheta;
            curPointLabels[3] = (j-1) + (i+1)*nTheta;

            face curFace(curPointLabels);
            curFace = curFace.reverseFace();
            faces[fI++] = curFace;

            // faces()[fI++] = face(curPointLabels);
        }

        labelList curPointLabels(4, -1);

        curPointLabels[0] = (nTheta-1) + i*nTheta;
        curPointLabels[1] = 0 + i*nTheta;
        curPointLabels[2] = 0 + (i+1)*nTheta;
        curPointLabels[3] = (nTheta-1) + (i+1)*nTheta;

        face curFace(curPointLabels);
        curFace = curFace.reverseFace();
        faces[fI++] = curFace;

        // faces()[fI++] = face(curPointLabels);
    }


    // Owners
    //Xfer<labelList> owners(labelList(nFaces, -1));
	labelList owners(nFaces, -1);

    fI = 0;
    for (label i=1; i<=nInternalFaces; i++)
    {
        owners[fI++] = i-1;
    }

    owners[fI++] = 0;

    owners[fI++] = nZ-1;

    for (label i=0; i<nZ; i++)
    {
        for (label j=1; j<nTheta; j++)
        {
            owners[fI++] = i;
        }

        owners[fI++] = i;
    }


    // Neighbours
    labelList neighbours(nInternalFaces, -1);

    fI = 0;
    for (label i=1; i<=nInternalFaces; i++)
    {
        neighbours[fI++] = i;
    }

    // Create mesh
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        std::move(points),
        std::move(faces),
        std::move(owners),
        std::move(neighbours)
    );

    List<polyPatch*> patches(3);

    // Inlet patch
    patches[0] =
        polyPatch::New
        (
            polyPatch::typeName,
            "left",
            1,
            nInternalFaces,
            0,
            mesh.boundaryMesh()
        ).ptr();

    // Outlet patch
    patches[1] =
        polyPatch::New
        (
            polyPatch::typeName,
            "right",
            1,
            nInternalFaces+1,
            1,
            mesh.boundaryMesh()
        ).ptr();

    // Pipe patch
    patches[2] =
        polyPatch::New
        (
            emptyPolyPatch::typeName,
            "beam",
            nTheta*nZ,
            nInternalFaces+2,
            2,
            mesh.boundaryMesh()
        ).ptr();

    mesh.addPatches(patches);

//    mesh.removeFiles();

    List<pointZone*> pz(0);
    List<faceZone*> fz(0);
    List<cellZone*> cz(1);

    labelList allCells(mesh.nCells(), -1);
    forAll(allCells, cellI)
    {
        allCells[cellI] = cellI;
    }

    cz[cz.size() - 1] =
        new cellZone
        (
            regionName,
            allCells,
            cz.size() - 1,
            mesh.cellZones()
        );

    mesh.addZones(pz, fz, cz);

    mesh.write();

    #include "printMeshSummary.H"

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
