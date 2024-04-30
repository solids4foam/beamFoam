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
    This class contains the mathematical formulae to:
    1. Construct a 3 x 3 rotation matrix from a 3-D rotation vector using the Rodrigues' formula (member functions 'rotationMatrix') 
    2. Construct a (pseudo-)vector from a specified rotation matrix using the Spurrier algorithm, wherein
       quaternions are calculated from the rotation matix and then using the quaternions, the corresponding
       3-D (pseudo)-rotation matrix is calculated (member functions 'pseudoVector')
    3. Interpolate rotation matrices between cell-centre and face-centre values of a
       beam cell-segment (member function 'interpolatedRotationMatrix')
    4. Interpolate rotational strains between cell-centre and face-centre values of a
       beam cell-segment (member function 'meanLineCurvature')
    5. Calculation of spatial derivative of the rotation matrix (member function 'rotationAngleDerivative')

Author
    Zeljko Tukovic, FSB Zagreb
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef pseudoVector_H
#define pseudoVector_H

#include "volFields.H"
#include "pointFields.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<volVectorField> pseudoVector(const volTensorField& rotationMatrix);
    
tmp<vectorField> pseudoVector(const tensorField& rotationMatrix);

vector pseudoVector(const tensor& R);

tmp<tensorField> rotationMatrix(const vectorField& rotatioAngle);
tmp<volTensorField> rotationMatrix(const volVectorField& rotatioAngle);
tmp<surfaceTensorField> rotationMatrix(const surfaceVectorField& rotatioAngle);

tensor rotationMatrix(const vector& rotationAngle);

tmp<surfaceTensorField> interpolateRotationMatrix
(
    const volTensorField& R
);

tmp<surfaceVectorField> meanLineCurvature(const volTensorField& R);
    
void interpolateRotationMatrix
(
    const volTensorField& R,
    surfaceTensorField& Rf,
    surfaceVectorField& snGradRotAngle
);

void interpolateRotationMatrix
(
    const beamModel& bm,
    const surfaceTensorField& Rf,
    volTensorField& R
);
    
tmp<surfaceVectorField> rotationAngleDerivative
(
    const volTensorField& rotatioMatrix
);
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
