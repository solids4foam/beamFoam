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
    Foam::morisonForce::interpolateWave

Description


SourceFiles
    interpolateWave.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef interpolateWave_H
#define interpolateWave_H

#include "morisonForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interpolateWave Declaration
\*---------------------------------------------------------------------------*/

class interpolateWave
:
    public morisonForce
{

protected:

    // Protected data
    scalarField z_;
    scalarField t_;
    vector invc_;

    vector velMean_;
    List<vectorField> velocity_;
    List<vectorField> acceleration_;

    // Protected member functions
    vector interpolateField
    (
        const List<vectorField>&,
        const scalar&,
        const scalar
    ) const;

public:

    //- Runtime type information
    TypeName("interpolateWave");

    // Constructors

        //- from components
        interpolateWave
        (
            const Time&,
            const fvMesh&,
            const dictionary&
        );


    // Destructor

        ~interpolateWave()
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
            notImplemented("Error!");
            return 1.0;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
