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
    stabilisedPenaltyFriction

Description
    Standard penalty method for normal contact law.

    Using a fixed penalty factor, force is applied to the contact faces based on
    the gap.

    The penalty method is effectively like attaching springs to any face that
    comes into contact.

SourceFiles
    stabilisedPenaltyFriction.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef stabilisedPenaltyFriction_H
#define stabilisedPenaltyFriction_H

#include "frictionContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class stabilisedPenaltyFriction Declaration
\*---------------------------------------------------------------------------*/

class stabilisedPenaltyFriction
:
    public frictionContactModel
{
    // Private data

        //- Contact law dictionary
        dictionary frictionContactModelDict_;

        //- Pentalty stiffness
        scalar epsilon_;

        //- Pentalty stiffness due to rotation
        scalar epsilonRot_;
    
        //- Friction coefficient
        scalar frictionCoeff_;

        //- Stabilistion region relative width
        scalar gCoeff_;

        //- Axial relaxation factor
        scalar axialRelaxationFactor_;
    
        //- transversal relaxation factor
        scalar transversalRelaxationFactor_;
    
    // Private Member Functions

        //- Disallow copy construct
        stabilisedPenaltyFriction(const stabilisedPenaltyFriction&);

        //- Disallow default bitwise assignment
        void operator=(const stabilisedPenaltyFriction&);

public:

    //- Runtime type information
    TypeName("stabilisedPenalty");


    // Constructors

        //- Construct from dictionary
        stabilisedPenaltyFriction
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~stabilisedPenaltyFriction()
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

        //- Return penalty stiffness
        scalar epsilon(const bool translation) const
        {
            if (translation)
            {
                return epsilon_;
            }
            else
            {
                return epsilonRot_;
            }
        }
    
        //- Return state of frictional contact
        virtual bool stick
        (
            const scalar gap,
            const scalar fcn = 0,
            const scalar offset = 0
        ) const;
    
        //- Return gap offset
        virtual scalar newGapOffset
        (
            const scalar gap,
            const scalar gapIncrement,
            const scalar fcn = 0,
            const scalar offset = 0,
            const bool translation = true
        ) const;

        virtual scalar axialRelaxationFactor() const
        {
            return axialRelaxationFactor_;
        }
    
        virtual scalar transversalRelaxationFactor() const
        {
            return transversalRelaxationFactor_;
        }    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
