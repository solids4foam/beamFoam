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
    tanhPenaltyFriction

\*---------------------------------------------------------------------------*/

#include "tanhPenaltyFriction.H"
#include "volFields.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

  defineTypeNameAndDebug(tanhPenaltyFriction, 0);
  addToRunTimeSelectionTable(frictionContactModel, tanhPenaltyFriction, dictionary);


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tanhPenaltyFriction::tanhPenaltyFriction
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
    frictionCoeff_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("frictionCoeff", 0)
    ),
    gCoeff_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("gCoeff", 0.2)
    )
{}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

scalar tanhPenaltyFriction::contactForce
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool translation
) const
{
    scalar a = epsilon_/(frictionCoeff_*fcn);

    scalar fc = frictionCoeff_*fcn*::tanh(a*gap);

    fc = epsilon_*gap;
    
    return fc;
}

scalar tanhPenaltyFriction::contactForceDerivative
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool translation
) const
{
    scalar a = epsilon_/(frictionCoeff_*fcn);

    // scalar dfc = epsilon_/sqr(::cosh(a*gap));

    scalar dfc = epsilon_;
    
    return dfc;
}

scalar tanhPenaltyFriction::contactForceDerivativeOverFcn
(
    const scalar gap,
    const scalar fcn,
    const scalar offset,
    const bool translation
) const
{
    scalar a = epsilon_/(frictionCoeff_*fcn);

    // scalar dfc =
    //     frictionCoeff_*::tanh(a*gap)
    //   - 1.0/(fcn*fcn*sqr(::cosh(a*gap)));

    scalar dfc = 0;
    
    return dfc;
}
    
scalar tanhPenaltyFriction::maxFrictionForce
(
    const scalar normalContactForce
) const
{
    return frictionCoeff_*normalContactForce;
}

scalar tanhPenaltyFriction::slipGapCorrection
(
    const scalar gap,
    const scalar forceCorrection
) const
{
    return forceCorrection/epsilon_;
}
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
