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

Author
    Zeljko Tukovic, FSB Zagreb

\*---------------------------------------------------------------------------*/

#include "beamHelperFunctions.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

label sharedFace(const cell& c0, const cell& c1)
{
    label sharedFaceIndex = -1;

    if (c0 == c1)
    {
        return sharedFaceIndex;
    }
    
    forAll(c0, f0I)
    {
        forAll(c1, f1I)
        {
            if (c0[f0I] == c1[f1I])
            {
                sharedFaceIndex = c0[f0I];
                break;
            }
        }

        if (sharedFaceIndex != -1)
        {
            break;
        }
    }

    return sharedFaceIndex;
}
  
label neiPatchIndex(const cell& c, const fvMesh& mesh)
{
    label patchIndex = -1;

    forAll(c, fI)
    {
        if (c[fI] >= mesh.nInternalFaces())
        {
            label patchI = mesh.boundaryMesh().whichPatch(c[fI]);

            if (!isA<emptyFvPatch>(mesh.boundary()[patchI]))
            {
                patchIndex = patchI;
                break;
            }
        }
    }

    return patchIndex;
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
