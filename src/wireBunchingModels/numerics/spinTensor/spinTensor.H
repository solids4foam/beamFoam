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
    This code is used for:
    1. Construction of a 3 x 3 skew-symmetric matrix, i.e., A.T() = - A,
       from a 3 x 1 vector (member function 'spinTensor').
    2. Construction of a 3 x 1 vector (also called as axial vector) from
       a skew-symmetric matrix (member function 'axialVector')

Author
    Zeljko Tukovic, FSB Zagreb,
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef spinTensor_H
#define spinTensor_H

#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<surfaceTensorField> spinTensor(const surfaceVectorField& axialVector);
tmp<volTensorField> spinTensor(const volVectorField& axialVector);
tmp<tensorField> spinTensor(const vectorField& axialVector);

tensor spinTensor(const vector& axialVector);

tmp<surfaceVectorField> axialVector(const surfaceTensorField& spinTensor);
tmp<volVectorField> axialVector(const volTensorField& spinTensor);
tmp<vectorField> axialVector(const tensorField& spintTensor);

vector axialVector(const tensor& spinTensor);
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
