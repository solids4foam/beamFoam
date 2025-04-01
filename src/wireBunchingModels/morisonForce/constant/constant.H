/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::morisonForce::constant

Description


SourceFiles
    constant.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "morisonForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace morison
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public morisonForce
{

protected:

    // Protected data
        vector Um_;
        scalar Tsoft_;

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- from components
        constant
        (
            const Time&,
            const fvMesh&,
            const dictionary&
        );


    // Destructor

        ~constant()
        {}


    // Member Functions
        virtual void printCoeffs();

        virtual dimensionedVector U
        (
            const point&,
            const scalar&
        ) const;

        virtual dimensionedVector acc
        (
            const point&,
            const scalar&
        ) const;

        scalar Uabsolute() const
        {
             return Foam::mag(Um_);
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace morison

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
