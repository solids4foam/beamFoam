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

#include "cylinderCellMarker.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Helper function to check if a point lies within a cylinder
bool isPointInCylinder(const point& p, const point& p1, const point& p2, double radius)
{
    vector d = p2 - p1;
    vector ap = p - p1;

    double dLength = mag(d);
    vector dUnit = d / dLength;

    // Project ap onto d to find the closest point along the cylinder axis
    double proj = (ap & dUnit);

    // Check if the projection is greater or smaller than the height of the cylinder
    if (proj < 0.0 || proj > dLength)
    {
        return false;
    }

    // Calculate the perpendicular distance from the point to the cylinder axis
    vector closestPoint = p1 + proj * dUnit;
    double distToAxis = mag(p - closestPoint);

    return distToAxis <= radius;
}

void markCellsInCylinders(
    const fvMesh& mesh,
    const List<point>& points,
    const scalar radius,
    volScalarField& cellMarker
)
{
    for (label i = 0; i < points.size() - 1; i++)
    {
        const point& p1 = points[i];
        const point& p2 = points[i + 1];
        
        forAll(mesh.C(),cellI)
        {
            if (isPointInCylinder(mesh.C()[cellI], p1, p2, radius))
            {
                cellMarker[cellI] = 1.0;
            }
        }
    }
}
} // End namspace Foam
