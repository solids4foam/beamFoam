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

Class
    tanhPenaltyFriction

Description
    Standard penalty method for normal contact law.

    Using a fixed penalty factor, force is applied to the contact faces based on
    the gap.

    The penalty method is effectively like attaching springs to any face that
    comes into contact.

SourceFiles
    tanhPenaltyFriction.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef tanhPenaltyFriction_H
#define tanhPenaltyFriction_H

#include "frictionContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tanhPenaltyFriction Declaration
\*---------------------------------------------------------------------------*/

class tanhPenaltyFriction
:
    public frictionContactModel
{
    // Private data

        //- Contact law dictionary
        dictionary frictionContactModelDict_;

        //- Pentalty stiffness
        scalar epsilon_;

        //- Friction coefficient
        scalar frictionCoeff_;

        //- Stabilistion region relative width
        scalar gCoeff_;
    
    // Private Member Functions

        //- Disallow copy construct
        tanhPenaltyFriction(const tanhPenaltyFriction&);

        //- Disallow default bitwise assignment
        void operator=(const tanhPenaltyFriction&);

public:

    //- Runtime type information
    TypeName("tanhPenalty");


    // Constructors

        //- Construct from dictionary
        tanhPenaltyFriction
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~tanhPenaltyFriction()
        {}


    // Member Functions

        //- Return normal contact force
        virtual scalar contactForce
        (
            const scalar gap,
            const scalar fcn,
            const scalar offset,
            const bool translation = true
        ) const;

        //- Return normal contact force deriative
        virtual scalar contactForceDerivative
        (
            const scalar gap,
            const scalar fcn,
            const scalar offset,
            const bool translation = true
        ) const;
  
        //- Return normal contact force deriative
        virtual scalar contactForceDerivativeOverFcn
        (
            const scalar gap,
            const scalar fcn,
            const scalar offset,
            const bool translation = true
        ) const;
    
        //- Return normal contact force
        virtual scalar maxFrictionForce
        (
            const scalar normalContactForce
        ) const;

        //- Return slip gap correction
        virtual scalar slipGapCorrection
        (
            const scalar gap,
            const scalar forceCorrection
        ) const;

        //- Return friction coefficient
        virtual scalar frictionCoeff() const
        {
            return frictionCoeff_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
