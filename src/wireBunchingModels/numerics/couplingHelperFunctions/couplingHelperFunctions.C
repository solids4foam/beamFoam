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

#include "couplingHelperFunctions.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // TODO : change the type with respect to output (velocities)
    tmp<vectorField> getFluidVelocity(
        const fvMesh& fluidMesh, 
        const fvMesh& mesh, 
        const vectorField& beamCellCenterCoord,
        labelList& seedCellIDs,
        const scalar groundZ,
        const meshSearch& searchEngine // flag for using octree or not
        )
    {
        tmp<vectorField> tresult(new vectorField(mesh.nCells(),vector::zero));
        vectorField& result = tresult.ref();
        const volVectorField& fluidVelocity = fluidMesh.lookupObject<volVectorField>("U");


        //Todo which proc, 
        // labelList procID(mesh.nCells(), -1);

        // todo flag for octree option
         


        forAll(mesh.C(),beamCellI)
        {
            // Info << "the beam.z() is  =" << beamCellCenterCoord[beamCellI].z() << endl;
            // Info << " gz =  " << groundZ << endl;
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
                    // Info << "seed for cell 0  "  << seedCellIDs[beamCellI] << endl;
                }
                else
                {
                    seedCellIDs[beamCellI] = seedCellIDs[beamCellI-1] ;
                    // Info << "seed for cell " << beamCellI << "= "  << seedCellIDs[beamCellI] << endl;
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
                // Todo procID[beamCellID] = Pstream::myProcID();
            }
            seedCellIDs[beamCellI] = fluidCellID;
        }
        // reduce(sumOp<vectorField>(), result);
        // reduce(maxOp<labelList>(), procID);
        return tresult;
        
    }

} // End namespace Foam


// ************************************************************************* //
