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
    frictionContactModel

Description
    Model for the behaviour of the friction contact

SourceFiles
    frictionContactModel.C
    newFrictionContactModel.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef frictionContactModel_H
#define frictionContactModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------* \
                         Class frictionContactModel Declaration
\*---------------------------------------------------------------------------*/

class frictionContactModel
{
    // Private data

        //- Name
        const word name_;
  
    // Private Member Functions

        //- Disallow copy construct
        frictionContactModel(const frictionContactModel&);

        //- Disallow default bitwise assignment
        void operator=(const frictionContactModel&);

public:

    //- Runtime type information
    TypeName("frictionContactModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            frictionContactModel,
            dictionary,
            (
                const word name,
                const dictionary& dict
            ),
            (
                name,
                dict
            )
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<frictionContactModel> New
        (
            const word& name,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        frictionContactModel
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~frictionContactModel()
        {}


    // Member Functions

        //- Return friction contact force
        virtual scalar contactForce
        (
            const scalar gap,
            const scalar fcn = 0,
            const scalar offset = 0,
            const bool translation = true
        ) const = 0;

        //- Return friction contact force deriative
        virtual scalar contactForceDerivative
        (
            const scalar gap,
            const scalar fcn = 0,
            const scalar offset = 0,
            const bool translation = true
        ) const = 0;
  
        //- Return friction contact force deriative
        virtual scalar contactForceDerivativeOverFcn
        (
            const scalar gap,
            const scalar fcn = 0,
            const scalar offset = 0,
            const bool translation = true
        ) const
        {
            return 0;
        }
    
        //- Return friction contact force
        virtual scalar maxFrictionForce
        (
            const scalar normalContactForce
        ) const = 0;
  
        //- Return friction contact force
        virtual scalar slipGapCorrection
        (
            const scalar gap,
            const scalar forceCorrection
        ) const = 0;
  
        //- Return friction coefficient
        virtual scalar frictionCoeff() const = 0;

        //- Return state of frictional contact
        virtual bool stick
        (
            const scalar gap,
            const scalar fcn = 0,
            const scalar offset = 0
        ) const
        {
            return true;
        }
    
        //- Return gap offset
        virtual scalar newGapOffset
        (
            const scalar gap,
            const scalar gapIncrement,
            const scalar fcn = 0,
            const scalar offset = 0,
            const bool translation = true
        ) const
        {
            return offset;
        }

        virtual scalar axialRelaxationFactor() const
        {
            return 1.0;
        }
    
        virtual scalar transversalRelaxationFactor() const
        {
            return 1.0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
