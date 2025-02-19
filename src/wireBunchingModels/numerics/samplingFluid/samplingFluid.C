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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    std::pair<tmp<volVectorField>, tmp<volScalarField>> getFluidVelocity(
        const fvMesh& fluidMesh,
        const fvMesh& mesh,
        const vectorField& beamCellCenterCoord,
        labelList& seedCellIDs,
        const scalar groundZ,
        const meshSearch& searchEngine // flag for using octree or not
        )
    {
        tmp<volVectorField> tresult
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
    tmp<volScalarField> tmarkerResults
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
        volVectorField& result = tresult.ref();
        volScalarField& markerResults = tmarkerResults.ref();

        const volVectorField& fluidVelocity = fluidMesh.lookupObject<volVectorField>("U");

        //Todo which proc,
        // labelList procID(mesh.nCells(), -1);



        forAll(mesh.C(),beamCellI)
        {
            if (beamCellCenterCoord[beamCellI].z() <= groundZ)
            {
                seedCellIDs[beamCellI] = -1;
                result[beamCellI] = vector::zero;
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
            if (seedCellIDs[beamCellI] != -1) // && procID[beamCellI] == Pstream::myProcID())
            {
                // Info << "using walk " << endl;
                fluidCellID = searchEngine.findCell(beamCellCenterCoord[beamCellI],seedCellIDs[beamCellI]);
            }
            if (fluidCellID == -1) // We should perform global search if no other procs found it.
            {
                // Info << "global " << endl;
                // perform a tree search
                fluidCellID = searchEngine.findCell(beamCellCenterCoord[beamCellI],-1,true);
            }
            if (fluidCellID == -1)
            {
                result[beamCellI] = vector::zero;

                //procID[beamCellID] = -1;
            }
            else
            {
                result[beamCellI] = fluidVelocity[fluidCellID];
                markerResults[fluidCellID] = 1.0;

                // Todo procID[beamCellID] = Pstream::myProcID();
            }
            seedCellIDs[beamCellI] = fluidCellID;
        }
        //- creating a link to beamProperties dictionary
        const dictionary& linkToBeamProperties = mesh.time().db().parent().lookupObject<dictionary>("beamProperties");

        List<point> points(mesh.nCells() + 1);
        label startPatchID = mesh.boundaryMesh().findPatchID
            (
                "left"
            );
        label endPatchID = mesh.boundaryMesh().findPatchID
            (
                "right"
            );

        forAll(points, i)
        {
            if (i == 0)
            {
                points[i] = mesh.boundaryMesh()[startPatchID].faceCentres()[0]
                  + mesh.lookupObject<surfaceVectorField>("refWf").boundaryField()[startPatchID][0]
                  + mesh.lookupObject<volVectorField>("W").boundaryField()[startPatchID][0];
            }
            else if (i == mesh.nCells())
            {
                points[i] = mesh.boundaryMesh()[endPatchID].faceCentres()[0]
                  + mesh.lookupObject<surfaceVectorField>("refWf").boundaryField()[endPatchID][0]
                  + mesh.lookupObject<volVectorField>("W").boundaryField()[endPatchID][0];
            }
            else
            {
                points[i] = mesh.Cf()[i-1]
                  + mesh.lookupObject<surfaceVectorField>("refWf")[i-1]
                  + mesh.lookupObject<volVectorField>("W")[i-1];
            }
        }
        //- access to beam radius
        dimensionedScalar radius = linkToBeamProperties.get<dimensionedScalar>("R");

        //- Calling the function to mark the cells in fluidMesh
        markCellsInCylinders(fluidMesh, points, radius.value(), markerResults);

        return std::make_pair(tresult, tmarkerResults);

    }

} // End namespace Foam


// ************************************************************************* //
