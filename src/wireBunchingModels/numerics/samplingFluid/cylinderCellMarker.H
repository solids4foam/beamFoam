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
    A function to mark the fluid cells that are lying inside a beam cells
Author
    Amirhossein Taran, UCD

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "point.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifndef cylinderCellMarker_H
#define cylinderCellMarker_H

namespace Foam
{
    
// Helper function to check if a point lies within a cylinder
bool isPointInCylinder(const point& p, const point& p1, const point& p2, double radius);

    
void markCellsInCylinders(
    const fvMesh& mesh,
    const List<point>& points,
    const scalar radius,
    //const List<point>& centers,
    volScalarField& cellMarker
);

} // End namespaceFoam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
