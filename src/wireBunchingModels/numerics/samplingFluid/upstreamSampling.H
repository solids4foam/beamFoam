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
    Helper function for Actuator line model's upstream sampling
Author
    Amirhossein Taran, UCD

\*---------------------------------------------------------------------------*/
#ifndef upstreamSampling_H
#define upstreamSampling_H

#include "fvMesh.H"
#include "vector.H"
#include "pointField.H"
#include "volFields.H"

namespace Foam
{

// Compute upstream direction perpendicular to tangent
inline vector upstreamDirPerpToTangent(const vector& vf, const vector& dRdScell)
{
    vector t = dRdScell/(mag(dRdScell) + VSMALL);
    vector vPerp = vf - (vf & t)*t;
    scalar m = mag(vPerp);
    if (m < VSMALL)
    {
        return -t;
    }
    return -vPerp/(m + VSMALL);
}

// Compute sampling points Ps (step 1 of ALM)
inline void computeUpstreamSamplingPoints(
    const pointField& beamCoords,
    const vectorField& beamUVec,
    const vectorField& fluidVelocity,
    const labelList& fluidCellIDs,
    const vectorField& dRdScell,
    const scalar sampleDist,
    pointField& Ps,
    vectorField& upstreamDir
)
{
    Ps.setSize(beamCoords.size(), vector::zero);
    upstreamDir.setSize(beamCoords.size(), vector::zero);

    forAll(beamCoords, i)
    {
        vector vf = vector::zero;
        if (fluidCellIDs[i] >= 0)
        {
            vf = fluidVelocity[fluidCellIDs[i]] ;
        }

        vector dir = upstreamDirPerpToTangent(vf, dRdScell[i]);
        upstreamDir[i] = dir;
        Ps[i] = beamCoords[i] + sampleDist * dir;
    }
}

} // End namespace Foam

#endif