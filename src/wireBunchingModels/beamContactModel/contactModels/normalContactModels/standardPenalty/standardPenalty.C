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
    standardPenalty

\*---------------------------------------------------------------------------*/

#include "standardPenalty.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(standardPenalty, 0);
  addToRunTimeSelectionTable(normalContactModel, standardPenalty, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

standardPenalty::standardPenalty
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
    )
{}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

scalar standardPenalty::contactForce
(
    const scalar gap,
    const scalar forceOffset,
    const scalar gapOffset
) const
{
    scalar fc = forceOffset - epsilon_*(gap - gapOffset);

    if (fc < 0)
    {
        fc = 0;
    }

    return fc;
}

scalar standardPenalty::contactForceDerivative
(
    const scalar gap,
    const scalar forceOffset,
    const scalar gapOffset
) const
{
    scalar dfc = -epsilon_;

    scalar fc = forceOffset - epsilon_*(gap -  gapOffset);

    if (fc < 0)
    {
        dfc = 0;
    }
    
    return dfc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
