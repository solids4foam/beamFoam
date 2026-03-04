/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "samplingFluid.H"
#include <tuple>
#include "surfaceMesh.H"
#include "cylinderCellMarker.H"
#include "surfaceFields.H"
#include "Pstream.H"
#include "PstreamReduceOps.H"
#include "FieldSumOp.H"
#include "upstreamSampling.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    std::tuple
    <
        tmp<volVectorField>,
        tmp<volScalarField>,
        tmp<volVectorField>,
        labelList
    >
    getFluidVelocity
    (
        const fvMesh& fluidMesh,
        const fvMesh& mesh,
        const vectorField& beamCellCenterCoord,
        labelList& seedCellIDs,
        const scalar groundZ,
        const bool groundContactActive,
        const vectorField& dRdScell,
        const scalar samplingRadius,
        const meshSearch& searchEngine // flag for using octree or not
    )
    {
    tmp<volVectorField> tresult // this is the fluid velocity stored on beam mesh
    (
        new volVectorField
        (
            IOobject
            (
                "fluidVelocity",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "zero", dimVelocity, vector::zero
            )
        )
    );
    tmp<volScalarField> tmarkerResults // this is the cellMarker stored on fluid mesh
    (
        new volScalarField
        (
            IOobject
            (
                "fluidCellMarker",
                fluidMesh.time().timeName(),
                fluidMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidMesh,
            dimensionedScalar(dimless, 0)
        )
    );
    tmp<volVectorField> tmarkerVelocity // this is the beam velocity stored on fluid mesh
    (
        new volVectorField
        (
            IOobject
            (
                "markerVelocity",
                fluidMesh.time().timeName(),
                fluidMesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidMesh,
            dimensionedVector
            (
                "zero", dimVelocity, vector::zero
            )
        )
    );
    volVectorField& result = tresult.ref();
    volScalarField& markerResults = tmarkerResults.ref();
    volVectorField& markerVelocity = tmarkerVelocity.ref();

    vectorField beamCoords = beamCellCenterCoord;
    vectorField beamTangents =  dRdScell;
    const volVectorField& fluidVelocity = fluidMesh.lookupObject<volVectorField>("U");

    //- this is the beam cell velocities stored on the beam mesh
    // const volVectorField& beamU = mesh.lookupObject<volVectorField>("U");

    vectorField beamUVec(beamCoords.size(), vector::zero);

    if (mesh.nCells() != 0)
    {
        forAll(beamUVec, i)
        {
            beamUVec[i] = mesh.lookupObject<volVectorField>("U")[i];   // copy raw vectors
        }
    }
    Pstream::broadcast(beamUVec);
    Pstream::broadcast(beamCoords);
    Pstream::broadcast(beamTangents);


    seedCellIDs.setSize(beamCoords.size(), -1);
    // resultVect is the fluid velocity on the beam mesh
    vectorField resultVect(beamCoords.size(), vector::zero);


    labelList fluidCellIDs(beamCoords.size(), -1);
    labelList sampleCellIDs(beamCoords.size(), -1);
    boundBox bb = fluidMesh.bounds();
    bb.inflate(0.01);
    forAll(beamCoords, beamCellI)
    {
        if (!bb.contains(beamCoords[beamCellI]))
        {
            // if its not in the bounding box, return zero for fluid's U
            resultVect[beamCellI] = vector::zero;
            seedCellIDs[beamCellI] = -1;
            // if its not within the bounding box, set the fluidCellID to -1
            fluidCellIDs[beamCellI] = -1;
            continue;
        }
        // this part only should activate when we have a mooring case!
        if (beamCoords[beamCellI].z() <= groundZ && groundContactActive)
        {
            seedCellIDs[beamCellI] = -1;
            resultVect[beamCellI] = vector::zero;
            continue;
        }
        if (seedCellIDs[beamCellI] == -1)
        {
            if (beamCellI == 0)
            {
                seedCellIDs[beamCellI] = seedCellIDs[beamCellI+1] ; // look for the seed in the next cell
            }
            else
            {
                seedCellIDs[beamCellI] = seedCellIDs[beamCellI-1] ; // look for the seed in the previous cell
            }
        }
        label fluidCellID = -1;
        if (seedCellIDs[beamCellI] != -1)
        {
            fluidCellID = searchEngine.findCell(beamCoords[beamCellI],seedCellIDs[beamCellI]);
        }

        if (fluidCellID == -1)
        {
            fluidCellID = searchEngine.findCell(beamCoords[beamCellI],-1,true);
        }
        // commented out since we are not getting the velocity here!
        if (fluidCellID != -1)
        {
            fluidCellIDs[beamCellI] = fluidCellID;

            markerVelocity[fluidCellID] = beamUVec[beamCellI];
        }
        // else
        // {
        // }

        // fluidcellID may or may not be -1
        // if -1, we could keep the old seed ...
        seedCellIDs[beamCellI] = fluidCellID;
    }
    // @todo check that each beam cell is found by only one procs.

    // Initializing radius on all ranks
    scalar radius(0.0);
    // ALM sampling
    if (Pstream::master())
    {
        const dictionary& linkToBeamProperties =
            mesh.time().db().parent().lookupObject<dictionary>
            (
                "beamProperties"
            );
        radius = linkToBeamProperties.get<dimensionedScalar>("R").value();
    }
    Pstream::broadcast(radius);

    pointField Ps;
    vectorField upstreamDir;

    const scalar sampleDist = samplingRadius*radius;

    vectorField vfAtBeam(beamCoords.size(), vector::zero);
    forAll(beamCoords, i)
    {
        if (fluidCellIDs[i] >= 0)
        {
            // local rank contributes the velocity from its own cell
            vector U_fluid = fluidVelocity[fluidCellIDs[i]];
            // We should check the beam UVec!
            vector U_beam = beamUVec[i];
            // get the relative velocity, and use it for getting the correct sampling directions
            vfAtBeam[i] = U_fluid - U_beam; // check here @todo
        }
    }
    // Reduce so that *all ranks* now have the SAME vfAtBeam
    // this assumes each beam cell maps to only one procs.
    reduce(vfAtBeam, FieldSumOp<vector>());


    //---- parallel code up until here
    if (Pstream::master())
    {
        computeUpstreamSamplingPoints(
            beamCoords,    // broadcast earlier
            vfAtBeam,      // fluid velocity at beam points, reduced across ranks
            beamTangents,  // broadcast earlier
            sampleDist,
            Ps,
            upstreamDir
        );
    }

    // Broadcast upstream sampling points and directions so all ranks use the same Ps / upstreamDir
    Pstream::broadcast(Ps);
    Pstream::broadcast(upstreamDir);


    forAll(Ps, pI)
    {
        label sampleCellID = -1;
        if (!fluidMesh.bounds().contains(Ps[pI]))
        {
            resultVect[pI] = vector::zero;
            continue;
        }
        // try to also add searching with seed
        sampleCellID = searchEngine.findCell(Ps[pI],-1,true);
        if (sampleCellID == -1)
        {
            const vector dirHat = upstreamDir[pI]/(mag(upstreamDir[pI]) + VSMALL);
            const point pPert = Ps[pI] + 1e-9*dirHat;   // small directional nudge
            sampleCellID = searchEngine.findCell(pPert, -1, true);
        }
        if (sampleCellID == -1)
        {
            resultVect[pI] = vector::zero;
        }
        else
        {
            resultVect[pI] = fluidVelocity[sampleCellID];
            markerResults[sampleCellID] = 1.0;
        }
    }


    reduce(resultVect, FieldSumOp<vector>());
    result.primitiveFieldRef() = resultVect;

    return std::make_tuple(std::move(tresult), std::move(tmarkerResults), std::move(tmarkerVelocity), std::move(fluidCellIDs));
    }
} // End namespace Foam


// ************************************************************************* //
