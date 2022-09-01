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
    stabilisedPenaltyFriction

\*---------------------------------------------------------------------------*/

#include "stabilisedPenaltyFriction.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(stabilisedPenaltyFriction, 0);
  addToRunTimeSelectionTable(frictionContactModel, stabilisedPenaltyFriction, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

stabilisedPenaltyFriction::stabilisedPenaltyFriction
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
    epsilonRot_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("epsilonRot", epsilon_)
    ),
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

scalar stabilisedPenaltyFriction::contactForce
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    scalar smallFcn = 0;

    // scalar gCoeff = 0;
    scalar gCoeff = gCoeff_;

    scalar FCN = fcn + smallFcn;

    scalar gStar = frictionCoeff_*FCN/epsilon(trans);
    scalar gBar = gCoeff*gStar;
    scalar gStart = gStar - gBar;
    scalar gEnd = gStar + gBar;

    scalar fc = frictionCoeff_*FCN;

    scalar relGap = gap - offset;

    scalar g = mag(relGap);

    if (g < gStart)
    {
        fc = epsilon(trans)*g;
    }
    else if (g < gEnd)
    {
        fc = (1.0-gCoeff)*frictionCoeff_*FCN
          + epsilon(trans)*(g - gStart)
          - (sqr(epsilon(trans))/(4*gCoeff*frictionCoeff_*FCN))
           *sqr(g - gStart);
    }

    // Set sign
    fc *= relGap/(mag(relGap) + SMALL);

    return fc;
}

scalar stabilisedPenaltyFriction::contactForceDerivative
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    scalar smallFcn = 0;

    scalar gCoeff = gCoeff_;
    // scalar gCoeff = 1.15*gCoeff_;

    scalar FCN = fcn + smallFcn;
    
    scalar gStar = frictionCoeff_*FCN/epsilon(trans);
    scalar gBar = gCoeff*gStar;
    scalar gStart = gStar - gBar;
    scalar gEnd = gStar + gBar;
    
    scalar dfc = 0;

    scalar relGap = gap - offset;

    scalar g = mag(relGap);
    scalar gCorrected = g;

    scalar dg = 0.2*gBar;
    // scalar curDq = (dg/(2*gBar))*(g - gStart);

    if (g > gStart && g <= gEnd)
    {
        scalar curDq = (dg/(2*gBar))*(g - gStart);
        gCorrected -= curDq;
    }
    else if (g > gEnd)
    {
        gCorrected -= dg;
    }

    // scalar gCorrected = g - curDq;

    if (gCorrected <= gStart)
    {
        dfc = epsilon(trans);
    }
    else if (gCorrected <= gEnd)
    {
        // dfc =
        //     epsilon(trans)
        //   - (epsilon(trans)/(gEndStab - 47gStart))*(g - gStart);
        dfc =
            epsilon(trans)
          - (sqr(epsilon(trans))
           /(2*gCoeff*frictionCoeff_*FCN))*(gCorrected - gStart);
    }

    return dfc;
}

scalar stabilisedPenaltyFriction::contactForceDerivativeOverFcn
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    scalar smallFcn = 0;
    
    // scalar gCoeff = 0;
    scalar gCoeff = gCoeff_;

    scalar FCN = fcn + smallFcn;

    scalar gStar = frictionCoeff_*FCN/epsilon(trans);
    scalar gBar = gCoeff*gStar;
    scalar gStart = gStar - gBar;
    scalar gEnd = gStar + gBar;

    scalar dfc = frictionCoeff_;

    scalar relGap = gap - offset;
    
    scalar g = mag(relGap);

    // if (g < gStart)
    // {
    //     dfc = 0;
    // }
    // else if (g < gEnd)
    // {
    //     dfc = (1.0 - gCoeff_)*frictionCoeff_
    //       - sqr(epsilon(trans))*sqr(g - gStart)
    //        /(4*gCoeff_*frictionCoeff_*sqr(fcn));
    // }

    if (g < gStart)
    {
        dfc = 0;
    }
    else if (g < gEnd)
    {
        dfc =
            (
                sqr(epsilon(trans))*(g - gStart)
               /(4*gCoeff*frictionCoeff_*FCN)
            )
           *(
               (g - gStart)/FCN
              + 2*(1.0 - gCoeff)
               *frictionCoeff_/epsilon(trans)
            );
    }

    dfc *= relGap/(mag(relGap) + SMALL);

    return dfc;
}
    
scalar stabilisedPenaltyFriction::maxFrictionForce
(
    const scalar normalContactForce
) const
{
    return frictionCoeff_*normalContactForce;
}

scalar stabilisedPenaltyFriction::slipGapCorrection
(
    const scalar gap,
    const scalar forceCorrection
) const
{
    return forceCorrection/epsilon_;
}

bool stabilisedPenaltyFriction::stick
(
    const scalar gap,
    const scalar fcn,
    const scalar offset
) const
{
    scalar smallFcn = 0;

    scalar gCoeff = gCoeff_;
    // scalar gCoeff = 0.5*gCoeff_;
    
    scalar FCN = fcn + smallFcn;

    // if (fcn < 1e-6)
    //     return true;

    scalar gStar = frictionCoeff_*FCN/epsilon_;
    scalar gBar = gCoeff*gStar;
    scalar gEnd = gStar + gBar;

    scalar relGap = gap - offset;

    scalar g = mag(relGap);

    if (g > gEnd)
    {
        return false;
    }

    return true;
}
    
scalar stabilisedPenaltyFriction::newGapOffset
(
    const scalar gap,
    const scalar gapIncrement,
    const scalar fcn,
    const scalar offset,
    const bool trans
) const
{
    scalar newOffset = offset;

    scalar smallFcn = 0;

    // scalar gCoeff = 0;
    scalar gCoeff = gCoeff_;

    scalar FCN = fcn + smallFcn;
    
    scalar gStar = frictionCoeff_*FCN/epsilon(trans);
    scalar gBar = gCoeff*gStar;
    scalar gStart = gStar - gBar;
    scalar gEnd = gStar + gBar;
    
    scalar relGap = gap - offset;
    
    scalar g = mag(relGap);

    // scalar relGapIncrement =
    //     gapIncrement/(relGap + SMALL);
    
    if (g > gEnd)
    {
        // scalar gTilda = gEnd - 0.1*gBar;
        // newOffset += (g - gTilda)*relGap/(mag(relGap) + SMALL);
        // newOffset += (g - gStar)*relGap/(mag(relGap) + SMALL);
        // newOffset += (g - gStart)*relGap/(mag(relGap) + SMALL);
        newOffset += (g - gEnd)*relGap/(mag(relGap) + SMALL);
    }
    // else if (relGapIncrement < 0)
    // {
    //     newOffset = gap - 0.5*(gap - offset);
    // }

    return newOffset;
}

    

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
