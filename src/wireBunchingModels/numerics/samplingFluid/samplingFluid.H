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

Description
    Helper function for two way coupled morrison equation
Author
    Amirhossein Taran, UCD
    Philip Cardiff, UCD

\*---------------------------------------------------------------------------*/

#ifndef samplingFluid_H
#define samplingFluid_H

#include "fvMesh.H"
#include "volFields.H"
#include "meshSearch.H"

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
            const meshSearch& searchEngine // flag for using octree
        );
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
