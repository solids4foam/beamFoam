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

    //- this field represents the beam cell velocity mapped on the fluid mesh
    volVectorField& markerVelocity = tmarkerVelocity.ref();

    vectorField beamCoords = beamCellCenterCoord;
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

    // commenting only for serial check
    Pstream::broadcast(beamCoords);


    seedCellIDs.setSize(beamCoords.size(), -1);
    vectorField resultVect(beamCoords.size(), vector::zero);

    globalIndex gI(fluidMesh.nCells());

    labelList fluidCellIDs(beamCoords.size(), -1);
    labelList sampleCellIDs(beamCoords.size(), -1);
    // markerResults.storePrevIter();
    forAll(beamCoords, beamCellI)
    {
        if (!fluidMesh.bounds().contains(beamCoords[beamCellI]))
        {
            resultVect[beamCellI] = vector::zero;
            seedCellIDs[beamCellI] = -1;
            // if its not within the bounding box, set the fluidCellID to -1!
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
                seedCellIDs[beamCellI] = seedCellIDs[beamCellI+1] ; // look for seed in the next cells
            }
            else
            {
                seedCellIDs[beamCellI] = seedCellIDs[beamCellI-1] ;
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
        // if (fluidCellID == -1)
        // {
        //     resultVect[beamCellI] = vector::zero;
        // }
        else
        {
            // fluidCellIDs[beamCellI] = fluidCellID;
            fluidCellIDs[beamCellI] = gI.toGlobal(fluidCellID);
            // resultVect[beamCellI] = fluidVelocity[fluidCellID];
            markerResults[fluidCellID] = 1.0;
            markerVelocity[fluidCellID] = beamUVec[beamCellI];
        }
        seedCellIDs[beamCellI] = fluidCellID;
    }
    // ALM sampling
    const dictionary& linkToBeamProperties =
        mesh.time().db().parent().lookupObject<dictionary>
        (
            "beamProperties"
        );
    dimensionedScalar radius = linkToBeamProperties.get<dimensionedScalar>("R");

    pointField Ps;
    vectorField upstreamDir;
    const scalar sampleDist = 2.0*radius.value();

    // Use seedCellIDs (local fluid cell indices) for accessing fluidVelocity
    computeUpstreamSamplingPoints(
        beamCoords,
        beamUVec,
        fluidVelocity,
        seedCellIDs,
        dRdScell,
        sampleDist,
        Ps,
        upstreamDir
    );
    label sampleCellID(-1);

    forAll(Ps, pI)
    {
        if (fluidCellIDs[pI] != -1)
        {
            sampleCellID = searchEngine.findCell(Ps[pI],fluidCellIDs[pI]);
        }
        if (sampleCellID == -1)
        {
            sampleCellID = searchEngine.findCell(Ps[pI],-1,true);
        }
        if (sampleCellID == -1)
        {
            resultVect[pI] = vector::zero;
        }
        else
        {
            resultVect[pI] = fluidVelocity[sampleCellID];
        }
    }

    // markerResults = alpha * markerResults + (1 - alpha)*markerResultsPrevIter;
    reduce(resultVect, FieldSumOp<vector>());

    forAll(fluidCellIDs, i)
    {
        reduce(fluidCellIDs[i], maxOp<label>());
    }
    result.primitiveFieldRef() = resultVect;

    // reduce(markerVelocityVect, FieldSumOp<vector>());
    // markerVelocity.primitiveFieldRef() = markerVelocityVect;
    return std::make_tuple(std::move(tresult), std::move(tmarkerResults), std::move(tmarkerVelocity), std::move(fluidCellIDs));
    }
} // End namespace Foam


// ************************************************************************* //
