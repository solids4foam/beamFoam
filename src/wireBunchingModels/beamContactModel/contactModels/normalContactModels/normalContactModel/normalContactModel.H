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
    normalContactModel

Description
    Model for the behaviour of the normal contact

SourceFiles
    normalContactModel.C
    newNormalContactModel.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef normalContactModel_H
#define normalContactModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------* \
                         Class normalContactModel Declaration
\*---------------------------------------------------------------------------*/

class normalContactModel
{
    // Private data

        //- Name
        const word name_;
  
    // Private Member Functions

        //- Disallow copy construct
        normalContactModel(const normalContactModel&);

        //- Disallow default bitwise assignment
        void operator=(const normalContactModel&);

public:

    //- Runtime type information
    TypeName("normalContactModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            normalContactModel,
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
        static autoPtr<normalContactModel> New
        (
            const word& name,
            const dictionary& dict
        );


    // Constructors

        //- Construct from dictionary
        normalContactModel
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~normalContactModel()
        {}


    // Member Functions

        //- Return normal contact force
        virtual scalar contactForce
        (
            const scalar gap,
            const scalar forceOffset = 0,
            const scalar gapOffset = 0
        ) const = 0;

        //- Return normal contact force deriative
        virtual scalar contactForceDerivative
        (
            const scalar gap,
            const scalar forceOffset = 0,
            const scalar gapOffset = 0
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
