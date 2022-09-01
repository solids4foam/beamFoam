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
    slipPenaltyFriction

\*---------------------------------------------------------------------------*/

#include "slipPenaltyFriction.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(slipPenaltyFriction, 0);
  addToRunTimeSelectionTable(frictionContactModel, slipPenaltyFriction, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

slipPenaltyFriction::slipPenaltyFriction
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
    frictionCoeff_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("frictionCoeff", 0)
    ),
    gCoeff_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("gCoeff", 0.2)
    ),
    axialRelaxationFactor_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("axialRelaxationFactor", 1.0)
    ),
    transversalRelaxationFactor_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("transversalRelaxationFactor", 1.0)
    )
{}//

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

scalar slipPenaltyFriction::contactForce
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    scalar fc = frictionCoeff_*fcn;

    scalar relGap = gap - offset;
    scalar g = mag(relGap);
    
    if (g<gCoeff_)
    {
        fc = 0;
    }
    
    // Set sign
    fc *= relGap/(mag(relGap) + SMALL);

    return fc;
}

scalar slipPenaltyFriction::contactForceDerivative
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    // scalar relGap = gap - offset;
    
    scalar dfc = 0;

    return dfc;
}

scalar slipPenaltyFriction::contactForceDerivativeOverFcn
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    scalar dfc = frictionCoeff_;

    scalar relGap = gap - offset;
    scalar g = mag(relGap);
    
    if (g<gCoeff_)
    {
        dfc = 0;
    }
    
    dfc *= relGap/(mag(relGap) + SMALL);

    return dfc;
}

scalar slipPenaltyFriction::maxFrictionForce
(
    const scalar normalContactForce
) const
{
    return frictionCoeff_*normalContactForce;
}

scalar slipPenaltyFriction::slipGapCorrection
(
    const scalar gap,
    const scalar forceCorrection
) const
{
    return 0;
}

bool slipPenaltyFriction::stick
(
    const scalar gap,
    const scalar fcn,
    const scalar offset
) const
{
    return false;
}
    
scalar slipPenaltyFriction::newGapOffset
(
    const scalar gap,
    const scalar gapIncrement,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    scalar newOffset = offset;

    scalar relGap = gap - offset;  
    scalar g = mag(relGap);

    if ((g-gCoeff_) > SMALL)
    {
        newOffset += (g-gCoeff_)*relGap/(mag(relGap) + SMALL);
    }

    return newOffset;
}

    

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
