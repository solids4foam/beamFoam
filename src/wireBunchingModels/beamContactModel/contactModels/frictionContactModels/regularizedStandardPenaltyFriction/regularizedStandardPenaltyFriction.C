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
    regularizedStandardPenaltyFriction

\*---------------------------------------------------------------------------*/

#include "regularizedStandardPenaltyFriction.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(regularizedStandardPenaltyFriction, 0);
  addToRunTimeSelectionTable(frictionContactModel, regularizedStandardPenaltyFriction, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

regularizedStandardPenaltyFriction::regularizedStandardPenaltyFriction
(
    const word& name,
    const dictionary& dict
)
:
    frictionContactModel
    (
        name,
        dict
    ),
    frictionContactModelDict_(dict.subDict(name + "FrictionModelDict")),
    epsilon_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("epsilon", 0)
    ),
    gBar_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("gBar", 0)
    ),
    frictionCoeff_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("frictionCoeff", 0)
    )
{}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

scalar regularizedStandardPenaltyFriction::contactForce
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool translation
) const
{
    scalar fc = 0;

    if ( gap > gBar_)
    {
        fc = epsilon_*(gap - gBar_/2);
    }
    else if (gap > 0)
    {
        fc = (epsilon_/(2*gBar_))*sqr(gap);
    }
  
    return fc;
}

scalar regularizedStandardPenaltyFriction::contactForceDerivative
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool translation
) const
{
    scalar dfc = 0;
  
    if ( gap > gBar_)
    {
        dfc = epsilon_;
    }
    else if (gap > 0)
    {
        dfc = (epsilon_/gBar_)*gap;
    }
    
    return dfc;
}

scalar regularizedStandardPenaltyFriction::maxFrictionForce
(
    const scalar normalContactForce
) const
{
    return frictionCoeff_*normalContactForce;
}

scalar regularizedStandardPenaltyFriction::slipGapCorrection
(
    const scalar gap,
    const scalar forceCorrection
) const
{
    scalar gapCorr = 0;

    scalar curForce = contactForce(gap, 0, 0);
    
    scalar newForce = curForce - forceCorrection;
        
    scalar newGap = 0;
    scalar fBar = epsilon_*gBar_/2;
        
    if (newForce >= fBar)
    {
        newGap = gBar_ + (newForce - fBar)/epsilon_;
    }
    else if (newForce > 0)
    {
        newGap = ::sqrt(2*gBar_*newForce/epsilon_);
    }

    gapCorr = gap - newGap;
    
    return gapCorr;
}
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
