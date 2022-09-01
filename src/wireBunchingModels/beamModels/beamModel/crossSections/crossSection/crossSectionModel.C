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

InClass
    crossSectionModel

\*---------------------------------------------------------------------------*/

#include "crossSectionModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(crossSectionModel, 0);
defineRunTimeSelectionTable(crossSectionModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
    
crossSectionModel::crossSectionModel
(
    const word& name,
    const dictionary& dict
)
:
    name_(name)
{}

    
void crossSectionModel::tensorToMatrix
(
    scalarRectangularMatrix& m,
    const tensor& t
) const
{
    m[0][0] = t.xx();
    m[0][1] = t.xy();
    m[0][2] = t.xz();

    m[1][0] = t.yx();
    m[1][1] = t.yy();
    m[1][2] = t.yz();

    m[2][0] = t.zx();
    m[2][1] = t.zy();
    m[2][2] = t.zz();
}


void crossSectionModel::tensorToMatrix
(
    scalarRectangularMatrix& m,
    const diagTensor& t
) const
{
    m[0][0] = t.xx();
    m[0][1] = 0;
    m[0][2] = 0;

    m[1][0] = 0;
    m[1][1] = t.yy();
    m[1][2] = 0;

    m[2][0] = 0;
    m[2][1] = 0;
    m[2][2] = t.zz();
}


void crossSectionModel::calcA
(
    scalarRectangularMatrix& A,
    const vector& p
) const
{
    A = 0*A;

    A[0][0] = 1;
    A[1][1] = 1;
    A[2][2] = 1;

    // A[0][4] =  p.y();
    // A[0][5] = -p.x();
    // A[1][3] = -p.y();
    // A[2][3] =  p.x();
    
    A[0][4] =  p.x();
    A[0][5] = -p.y();
    A[1][3] = -p.x();
    A[2][3] =  p.y();
}

scalar crossSectionModel::polarAngle(const scalar x, const scalar y) const
{
    scalar pAng = ::atan2(x, y);

    if (pAng < 0)
    {
        pAng = 2*M_PI + pAng;
    }

    return pAng;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
