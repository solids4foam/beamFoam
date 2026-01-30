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
// inline vector upstreamDirPerpToTangent(const vector& vf, const vector& beamTangents)
// {
//     vector t = beamTangents/(mag(beamTangents) + VSMALL);
//     vector vPerp = vf - (vf & t)*t;

<<<<<<< HEAD
    // extra checks to enforce upstream always to be in -X direction
    vector samplingDir = -vPerp;
    // // i just comment to check
    // if (samplingDir.x() > 0)
    // {
    //     vector globalUpstream(-1, 0, 0); 
    //     samplingDir = globalUpstream - (globalUpstream & t)*t;
    // }
=======
//     // extra checks to enforce upstream always to be in -X direction
//     vector samplingDir = -vPerp;
//     // if (samplingDir.x() > 0)
//     // {
//     //     vector globalUpstream(-1, 0, 0); 
//     //     samplingDir = globalUpstream - (globalUpstream & t)*t;
//     // }
>>>>>>> 95228d7 ( Refactor force calculation formulation and improved upstream direction handling in samplingFluid)
    
//     //    scalar m = mag(vPerp);
//     scalar m = mag(samplingDir);
//     if (m < VSMALL)
//     {
//         return -t;
//     }
//     //    return -vPerp/(m + VSMALL);
//     return samplingDir/(m + VSMALL);
// }
inline vector makePerpUnit(const vector& t)
{
    // Pick a reference axis that is NOT parallel to t
    const vector ex(1,0,0), ey(0,1,0), ez(0,0,1);

    vector a = ex;
    if (mag(t & ex) > 0.9) a = ey;      // if nearly parallel to x, use y
    if (mag(t & ey) > 0.9) a = ez;      // if nearly parallel to y, use z

    // Project a to be perpendicular to t
    vector p = a - (a & t)*t;

    // If still tiny (pathological), just use another axis
    if (mag(p) < VSMALL)
    {
        a = ez;
        p = a - (a & t)*t;
    }

    return p/(mag(p) + VSMALL);
}

inline vector upstreamDirPerpToTangent(const vector& vf, const vector& beamTangent)
{
    const vector t = beamTangent/(mag(beamTangent) + VSMALL);

    // Component perpendicular to tangent
    const vector vPerp = vf - (vf & t)*t;

    // Upstream = opposite of local perpendicular flow component
    vector dir = -vPerp;

    // If flow gives no useful direction, step away from the beam in a stable perpendicular direction
    if (mag(dir) < 1e-12)
    {
        dir = makePerpUnit(t);
    }

    return dir/(mag(dir) + VSMALL);
}

// Compute sampling points Ps (step 1 of ALM)
inline void computeUpstreamSamplingPoints(
    const pointField& beamCoords,
    const vectorField& vfAtBeam,
    const vectorField& beamTangents,
    const scalar sampleDist,
    pointField& Ps,
    vectorField& upstreamDir
)
{
    Ps.setSize(beamCoords.size(), vector::zero);
    upstreamDir.setSize(beamCoords.size(), vector::zero);

    forAll(beamCoords, i)
    {
        vector vf = vfAtBeam[i];

        vector dir = upstreamDirPerpToTangent(vf, beamTangents[i]);
        upstreamDir[i] = dir;
        Ps[i] = beamCoords[i] + sampleDist * dir;
    }
}

} // End namespace Foam

#endif
