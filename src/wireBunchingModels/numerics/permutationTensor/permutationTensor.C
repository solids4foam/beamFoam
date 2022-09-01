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

#include "permutationTensor.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<surfaceTensorField> mulByPermutationTensor
(
    const surfaceVectorField& svf,
    const bool right
)
{
    tmp<surfaceTensorField> tresult
    (
        new surfaceTensorField
        (
            IOobject
            (
                "permutationTensor",
                svf.time().timeName(),
                svf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            svf.mesh(),
            dimensionedTensor("0", svf.dimensions(), tensor::zero)
        )
    );

    surfaceTensorField& result = tresult();
    tensorField& resultI = result.internalField();

    const vectorField& svfI = svf.internalField();

    forAll(resultI, faceI)
    {
        resultI[faceI] = mulByPermutationTensor(svfI[faceI], right);
    }

    forAll(result.boundaryField(), patchI)
    {
        if
        (
            svf.boundaryField()[patchI].type()
         != emptyFvPatchField<vector>::typeName
        )
        {
            tensorField& pResult = result.boundaryField()[patchI];
            const vectorField& pSvf = svf.boundaryField()[patchI];

            forAll(pResult, faceI)
            {
                pResult[faceI] = mulByPermutationTensor(pSvf[faceI], right);
            }
        }
    }

    // result.correctBoundaryConditions();

    return tresult;
}

tensor mulByPermutationTensor(const vector& v, const bool right)
{
    tensor result = tensor::zero;

    tensor E0
    (
        0,  0, 0,
        0,  0, 1,
        0, -1, 0
    );

    tensor E1
    (
        0, 0, -1,
        0, 0,  0,
        1, 0,  0
    );

    tensor E2
    (
        0, 1, 0,
       -1, 0, 0,
        0, 0, 0
    );

    result.xx() = -(E0.xx()*v.x() + E0.xy()*v.y() + E0.xz()*v.z()); 
    result.xy() = -(E0.yx()*v.x() + E0.yy()*v.y() + E0.yz()*v.z()); 
    result.xz() = -(E0.zx()*v.x() + E0.zy()*v.y() + E0.zz()*v.z());
        
    result.yx() = -(E1.xx()*v.x() + E1.xy()*v.y() + E1.xz()*v.z()); 
    result.yy() = -(E1.yx()*v.x() + E1.yy()*v.y() + E1.yz()*v.z()); 
    result.yz() = -(E1.zx()*v.x() + E1.zy()*v.y() + E1.zz()*v.z());

    result.zx() = -(E2.xx()*v.x() + E2.xy()*v.y() + E2.xz()*v.z()); 
    result.zy() = -(E2.yx()*v.x() + E2.yy()*v.y() + E2.yz()*v.z()); 
    result.zz() = -(E2.zx()*v.x() + E2.zy()*v.y() + E2.zz()*v.z());

    if (right)
    {
        result *= -1;
    }

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
