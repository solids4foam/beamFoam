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
    circle

\*---------------------------------------------------------------------------*/

#include "circle.H"

#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(circle, 0);
  addToRunTimeSelectionTable(crossSectionModel, circle, dictionary);

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

circle::circle
(
    const word& name,
    const dictionary& dict
)
:
    crossSectionModel
    (
        name,
        dict
    ),
    crossSectionModelDict_(dict.subDict(name + "CrossSectionModelDict")),
    radius_
    (
        readScalar(crossSectionModelDict_.lookup("radius"))
    ),
    scalingArea_
    (
        crossSectionModelDict_.lookupOrDefault<scalar>("scalingArea", 1.0)
    ),
    scalingSecondMomentArea_
    (
        crossSectionModelDict_
        .lookupOrDefault<scalar>("scalingSecondMomentArea", 1.0)
    ),
    points_(0, vector::zero)
{}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

// tmp<vectorField> circle::greenLagrangianStrain
// (
//     const vector& Gamma,
//     const vector& K
// ) const
// {
//     tmp<vectorField> tE
//     (
//         new vectorField(points_.size(), vector::zero)
//     );

//     return tE;
// }


} // End namespace Foam

// ************************************************************************* //
