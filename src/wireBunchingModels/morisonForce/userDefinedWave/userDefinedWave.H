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
    Foam::morisonForce::userDefinedWave

Description


SourceFiles
    userDefinedWave.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef userDefinedWave_H
#define userDefinedWave_H

#include "morisonForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class userDefinedWave Declaration
\*---------------------------------------------------------------------------*/

class userDefinedWave
:
    public morisonForce
{

protected:

    // Protected data
    label N_;

    scalarField uSin_, uCos_, wSin_, wCos_;

    scalar u0_, du0dz_;
    scalar h_;
    scalar period_;
    scalar omega_;

    vector k_;
    scalar K_;

    scalar Tsoft_;

    Switch includeVertical_;
    Switch includeWaveNumber_;

    scalar spaceFact_;

    // Protected member functions
    void writeInputFile();

    void setCoefficients();

public:

    //- Runtime type information
    TypeName("userDefinedWave");

    // Constructors

        //- from components
        userDefinedWave
        (
            const Time&,
            const fvMesh&,
            const dictionary&
        );


    // Destructor

        ~userDefinedWave()
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
            Info << "Uabsolute needs implementing" << endl;
            return 1.0;
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
