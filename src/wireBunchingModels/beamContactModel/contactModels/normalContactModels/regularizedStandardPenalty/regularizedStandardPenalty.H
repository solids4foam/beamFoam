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
    regularizedStandardPenalty

Description
    Standard penalty method for normal contact law.

    Using a fixed penalty factor, force is applied to the contact faces based on
    the gap.

    The penalty method is effectively like attaching springs to any face that
    comes into contact.

SourceFiles
    regularizedStandardPenalty.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef regularizedStandardPenalty_H
#define regularizedStandardPenalty_H

#include "normalContactModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class regularizedStandardPenalty Declaration
\*---------------------------------------------------------------------------*/

class regularizedStandardPenalty
:
    public normalContactModel
{
    // Private data

        //- Contact law dictionary
        dictionary normalContactModelDict_;

        //- Pentalty stiffness
        scalar epsilon_;

        //- Regularization gap
        scalar gBar_;
  
    // Private Member Functions

        //- Disallow copy construct
        regularizedStandardPenalty(const regularizedStandardPenalty&);

        //- Disallow default bitwise assignment
        void operator=(const regularizedStandardPenalty&);

public:

    //- Runtime type information
    TypeName("regularizedStandardPenalty");


    // Constructors

        //- Construct from dictionary
        regularizedStandardPenalty
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~regularizedStandardPenalty()
        {}


    // Member Functions

        //- Return normal contact force
        virtual scalar contactForce
        (
            const scalar gap,
            const scalar forceOffset,
            const scalar gapOffset
        ) const;

        //- Return normal contact force deriative
        virtual scalar contactForceDerivative
        (
            const scalar gap,
            const scalar forceOffset,
            const scalar gapOffset
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
