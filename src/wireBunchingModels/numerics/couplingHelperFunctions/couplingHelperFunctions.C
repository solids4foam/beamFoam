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
    void amazingFunc(const fvMesh& fluidMesh, const fvMesh& mesh, const vectorField beamTotalDisp)
    {
        const volVectorField& fluidVelocity = fluidMesh.lookupObject<volVectorField>("U");
        meshSearch searchEngine(fluidMesh);
        DynamicList<label> closestFluidCellID(0);
        forAll(mesh.C(),beamCellI)
        {
            Info << "the location is  =" << beamTotalDisp[beamCellI] << endl;

            closestFluidCellID.append(searchEngine.findCell(beamTotalDisp[beamCellI],0,false));  
            // closestFluidCellID.append(fluidMesh.findCell(beamTotalDisp[beamCellI]));  
        }
    }

} // End namespace Foam


// ************************************************************************* //
