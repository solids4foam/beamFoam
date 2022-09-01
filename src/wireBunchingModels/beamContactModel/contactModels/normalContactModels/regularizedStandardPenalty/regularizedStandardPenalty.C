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
    regularizedStandardPenalty

\*---------------------------------------------------------------------------*/

#include "regularizedStandardPenalty.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(regularizedStandardPenalty, 0);
  addToRunTimeSelectionTable(normalContactModel, regularizedStandardPenalty, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regularizedStandardPenalty::regularizedStandardPenalty
(
    const word& name,
    const dictionary& dict
)
:
    normalContactModel
    (
        name,
        dict
    ),
    normalContactModelDict_(dict.subDict(name + "NormalModelDict")),
    epsilon_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("epsilon", 1.0)
    ),
    gBar_
    (
        normalContactModelDict_.lookupOrDefault<scalar>("gBar", 0)
    )
{}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

scalar regularizedStandardPenalty::contactForce
(
    const scalar gap,
    const scalar forceOffset,
    const scalar gapOffset
) const
{
    scalar fc = 0;

    scalar fBar = epsilon_*gBar_/2;

    scalar g = gap - gapOffset;
    
    if (g <= 0)
    {
        fc = fBar - epsilon_*g;
    }
    else if (g <= gBar_)
    {
        fc = (epsilon_*gBar_ - fBar)*sqr(g/gBar_)
          - epsilon_*g + fBar;
    }

    return fc;
}

scalar regularizedStandardPenalty::contactForceDerivative
(
    const scalar gap,
    const scalar forceOffset,
    const scalar gapOffset
) const
{
    scalar g = gap - gapOffset;
    
    scalar fBar = epsilon_*gBar_/2;

    scalar dfc = 0;

    if (g <= 0)
    {
        dfc = -epsilon_;
    }
    else if (g <= gBar_)
    {
        dfc = 2*g*(epsilon_*gBar_ - fBar)/sqr(gBar_)
          - epsilon_;
    }

    return dfc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
