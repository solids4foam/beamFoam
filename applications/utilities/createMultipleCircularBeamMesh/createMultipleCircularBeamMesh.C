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

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
    argList::validOptions.insert("decomposePatches", "");

#   include "setRootCase.H"
#   include "createTime.H"

    word regionName;
    fileName polyMeshDir;

    if (args.optionFound("region"))
    {
        // constant/<region>/polyMesh/blockMeshDict
        regionName  = args.option("region");
        polyMeshDir = regionName/polyMesh::meshSubDir;

        Info<< nl << "Generating mesh for region " << regionName << endl;
    }
    else
    {
        // constant/polyMesh/blockMeshDict
        regionName  = polyMesh::defaultRegion;
        polyMeshDir = polyMesh::meshSubDir;
    }

    bool decomposePatches = args.optionFound("decomposePatches");

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
            IOobject::NO_WRITE
        )
    );

    int nBeams = -1;
    scalarField R;
    scalarField L;
    labelList nZ;
    if (beamProperties.found("nBeams"))
    {
        nBeams = readInt(beamProperties.lookup("nBeams"));
        R.setSize
        (
            nBeams,
            dimensionedScalar(beamProperties.lookup("R")).value()
        );
        L.setSize
        (
            nBeams,
            dimensionedScalar(beamProperties.lookup("L")).value()
        );
        nZ.setSize
        (
            nBeams,
            readInt(beamProperties.lookup("nSegments"))
        );
    }
    else
    {
        R = scalarField(beamProperties.lookup("R"));
        L = scalarField(beamProperties.lookup("L"));
        nZ = labelList(beamProperties.lookup("nSegments"));
        if (R.size() != L.size())
        {
            FatalErrorIn(args.executable())
                << "Sizes of R and L fields are not equal."
                << abort(FatalError);
        }
        
        if (R.size() != nZ.size())
        {
            FatalErrorIn(args.executable())
                << "Sizes of R and nZ fields are not equal."
                << abort(FatalError);
        }
        
        nBeams = R.size();
    }

    if (nBeams < 1)
    {
        FatalErrorIn(args.executable())
            << "Number of beams is not defined."
            << abort(FatalError);
    }
    
    // scalar R = dimensionedScalar(beamProperties.lookup("R")).value();
    // scalar L = dimensionedScalar(beamProperties.lookup("L")).value();
    // scalar R = 0.1;
    // scalar L = 20;

    // label nZ(readInt(beamProperties.lookup("nSegments")));
    // label nZ = 50;
    label nTheta = 32;

    // const int nBeams
    // (
    //     beamProperties.lookupOrDefault<int>("nBeams", 1)
    // );

    labelList nOneBeamCells = nZ;

    labelList cellStart(nBeams, -1);
    cellStart[0] = 0;
    for (label i=1; i<nBeams; i++)
    {
        cellStart[i] = cellStart[i-1] + nOneBeamCells[i-1];
    }
    
    labelList nOneBeamPoints = nTheta*(nZ + 1);
    label nPoints = sum(nOneBeamPoints); //*nBeams;

    labelList pointStart(nBeams, -1);
    pointStart[0] = 0;
    for (label i=1; i<nBeams; i++)
    {
        pointStart[i] = pointStart[i-1] + nOneBeamPoints[i-1];
    }
    
    labelList nOneBeamFaces = ((nZ+1) + nZ*nTheta);
    label nFaces = sum(nOneBeamFaces); //*nBeams;

    labelList nOneBeamInternalFaces = (nZ-1);
    label nInternalFaces = sum(nOneBeamInternalFaces); //*nBeams;

    // Points
    Xfer<vectorField> points(vectorField(nPoints, vector::zero));

    scalar dTheta = -2*M_PI/nTheta;

    vectorField circumferencePoints(nTheta, vector::zero);
    forAll(circumferencePoints, pointI)
    {
        circumferencePoints[pointI].x() = ::cos(pointI*dTheta);
        circumferencePoints[pointI].y() = ::sin(pointI*dTheta);
        // circumferencePoints[pointI].x() = R*::cos(pointI*dTheta);
        // circumferencePoints[pointI].y() = R*::sin(pointI*dTheta);
        circumferencePoints[pointI].z() = 0;
    }

    label pI = 0;
    for (label k=0; k<nBeams; k++)
    {
        scalar dZ = L[k]/nZ[k];
        for (label i=0; i<=nZ[k]; i++)
        {
            for (label j=0; j<nTheta; j++)
            {
                vector curPoint = R[k]*circumferencePoints[j];
                curPoint.z() = i*dZ;

                points()[pI++] = curPoint;
            } 
        }
    }

    // Faces
    Xfer<faceList> faces(faceList(nFaces, face()));

    //-Internal faces

    label fI = 0;
    for (label k=0; k<nBeams; k++)
    {
        for (label i=0; i<=nZ[k]; i++)
        {
            labelList curPointLabels(nTheta, -1);
            for (label j=0; j<nTheta; j++)
            {
                curPointLabels[j] = j + i*nTheta + pointStart[k];
                // curPointLabels[j] = j + i*nTheta + k*nOneBeamPoints;
            }

            if (i>0 && i<nZ[k])
            {
                face curFace(curPointLabels);
                curFace = curFace.reverseFace();    
                faces()[fI++] = curFace;
            }
        }
    }

    //- Inlet faces

    for (label k=0; k<nBeams; k++)
    {
        labelList curPointLabels(nTheta, -1);
        for (label j=0; j<nTheta; j++)
        {
            curPointLabels[j] = j + pointStart[k]; //k*nOneBeamPoints;
        }
        
        face curFace(curPointLabels);

        faces()[fI++] = curFace;
    }
    
    //- Outlet faces

    for (label k=0; k<nBeams; k++)
    {
        labelList curPointLabels = labelList(nTheta, -1);
        for (label j=0; j<nTheta; j++)
        {
            curPointLabels[j] = j + nZ[k]*nTheta + pointStart[k]; //k*nOneBeamPoints;
        }

        face curFace = face(curPointLabels);
        curFace = curFace.reverseFace();    
        faces()[fI++] = curFace;
    }

    //- beam faces

    for (label k=0; k<nBeams; k++)
    {
        for (label i=0; i<nZ[k]; i++)
        {
            for (label j=1; j<nTheta; j++)
            {
                labelList curPointLabels(4, -1);

                curPointLabels[0] = (j-1) + i*nTheta;
                curPointLabels[1] = j + i*nTheta;
                curPointLabels[2] = j + (i+1)*nTheta;
                curPointLabels[3] = (j-1) + (i+1)*nTheta;

                forAll(curPointLabels, pI)
                {
                    curPointLabels[pI] += pointStart[k]; //k*nOneBeamPoints;
                }


                face curFace(curPointLabels);
                curFace = curFace.reverseFace();    
                faces()[fI++] = curFace;
            }

            labelList curPointLabels(4, -1);

            curPointLabels[0] = (nTheta-1) + i*nTheta;
            curPointLabels[1] = 0 + i*nTheta;
            curPointLabels[2] = 0 + (i+1)*nTheta;
            curPointLabels[3] = (nTheta-1) + (i+1)*nTheta;

            forAll(curPointLabels, pI)
            {
                curPointLabels[pI] += pointStart[k]; //k*nOneBeamPoints;
            }
      
            face curFace(curPointLabels);
            curFace = curFace.reverseFace();    
            faces()[fI++] = curFace;
        }
    }

    // Owners
    Xfer<labelList> owners(labelList(nFaces, -1));

    fI = 0;
    for (label k=0; k<nBeams; k++)
    {
        for (label i=0; i<=nZ[k]; i++)
        {
            if (i>0 && i<nZ[k])
            {
                owners()[fI++] = (i-1) + cellStart[k]; //k*nOneBeamCells;
            }
        }
    }

    // Inlet faces
    for (label k=0; k<nBeams; k++)
    {
        owners()[fI++] = cellStart[k]; //k*nOneBeamCells;
    }
    
    // Outlet faces
    for (label k=0; k<nBeams; k++)
    {
        owners()[fI++] = (nZ[k]-1) + cellStart[k]; //k*nOneBeamCells;
    }
    
    // Beam faces
    for (label k=0; k<nBeams; k++)
    {
        for (label i=0; i<nZ[k]; i++)
        {
            for (label j=1; j<nTheta; j++)
            {
                owners()[fI++] = i + cellStart[k]; //k*nOneBeamCells;            
            }

            owners()[fI++] = i + cellStart[k]; //k*nOneBeamCells;            
        }
    }
    
    // Neighbours
    Xfer<labelList> neighbours(labelList(nInternalFaces, -1));
    
    fI = 0;
    for (label k=0; k<nBeams; k++)
    {
        for (label i=1; i<=(nZ[k]-1); i++)
        {
            neighbours()[fI++] = i + cellStart[k]; //k*nOneBeamCells;
        }
    }
    
    // Create mesh
    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime
        ),
        points,
        faces,
        owners,
        neighbours
    );

    if ((nBeams > 1) && decomposePatches)
    {
        List<polyPatch*> patches(nBeams + nBeams + 1);

        label pI = 0;
        label start = nInternalFaces;
        
        // Inlet patches
        
        for (label k=0; k<nBeams; k++)
        {
            OStringStream PatchName;
            PatchName() << "left_" << k;

            patches[pI] =
                polyPatch::New
                (
                    polyPatch::typeName,
                    word(PatchName.str()),
                    1,
                    start,
                    pI,
                    mesh.boundaryMesh()
                ).ptr();
            pI++;
            start += 1;
        }

        // Outlet patches
        for (label k=0; k<nBeams; k++)
        {
            OStringStream PatchName;
            PatchName() << "right_" << k;
        
            patches[pI] = 
                polyPatch::New
                (
                    polyPatch::typeName,
                    word(PatchName.str()),
                    1,
                    start,
                    pI,
                    mesh.boundaryMesh()
                ).ptr();
            pI++;
            start += 1;
        }
        
        // Beam patch
        patches[pI] = 
            polyPatch::New
            (
                emptyPolyPatch::typeName,
                "beam",
                nTheta*sum(nZ),
                // (nTheta*nZ)*nBeams,
                start,
                pI,
                mesh.boundaryMesh()
            ).ptr();

        mesh.addPatches(patches);
    }
    else
    {
        List<polyPatch*> patches(3);

        // Inlet patch
        patches[0] = 
            polyPatch::New
            (
                polyPatch::typeName,
                "left",
                nBeams,
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
                nBeams,
                nInternalFaces + nBeams,
                1,
                mesh.boundaryMesh()
            ).ptr();

        // Beam patch
        patches[2] = 
            polyPatch::New
            (
                emptyPolyPatch::typeName,
                "beam",
                nTheta*sum(nZ),
                // (nTheta*nZ)*nBeams,
                nInternalFaces + nBeams + nBeams,
                2,
                mesh.boundaryMesh()
            ).ptr();

        mesh.addPatches(patches);
    }
    
    mesh.removeZones();

    List<pointZone*> pz(nBeams);
    List<faceZone*> fz(0);
    List<cellZone*> cz(nBeams);

    for (label k=0; k<nBeams; k++)
    {
        labelList curBeamCells(nOneBeamCells[k], -1);
        labelList curBeamPoints(nOneBeamPoints[k], -1);
        forAll(curBeamCells, cellI)
        {
            curBeamCells[cellI] = cellI + cellStart[k]; //k*nOneBeamCells;
        }
        forAll(curBeamPoints, pointI)
        {
            curBeamPoints[pointI] = pointI + pointStart[k]; //k*nOneBeamCells;
        }

        OStringStream RegionName;
        RegionName() << "beam_" << k;

        pz[k] =
            new pointZone
            (
                word(RegionName.str()),
                curBeamPoints,
                k,
                mesh.pointZones()
            );
        
        cz[k] =
            new cellZone
            (
                word(RegionName.str()),
                curBeamCells,
                k,
                mesh.cellZones()
            );
    }

    mesh.addZones(pz, fz, cz);
    
    mesh.write();

    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
