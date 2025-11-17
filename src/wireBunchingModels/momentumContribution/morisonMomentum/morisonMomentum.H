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
    Foam::momentumContribution::morisonMomentum

Description
    Derived class to add the Morison force contributions to the
    linear and angular momentum equations of the beam solver

SourceFiles
    morisonMomentum.C

Author
    Seevani Balii, UCD.
    Philip Cardiff, UCD.


\*---------------------------------------------------------------------------*/

#ifndef morisonMomentum_H
#define morisonMomentum_H

#include "momentumContribution.H"
#include "HermiteSpline.H"
// #include "beamModel.H"
#include "coupledTotalLagNewtonRaphsonBeam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// namespace momentumContributions
// {

/*---------------------------------------------------------------------------*\
                           Class morisonMomentum Declaration
\*---------------------------------------------------------------------------*/

class morisonMomentum
:
    public momentumContribution
{

private:

    //- Normal drag coefficient for Morison Force
    scalar Cdn_;

    //- Tangential drag coefficient for Morison Force
    scalar Cdt_;

public:

    //- Runtime type information
    TypeName("morisonMomentum");

    // Constructors

        //- Contruct from components
        morisonMomentum
        (
            const Time&
        );


    // Destructor

        ~morisonMomentum()
        {}


    // Member Functions

        //- Return the contribution to the block diagonal of the Jacobian matrix
        //  Each block should be a 6x6 matrix
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const beamModel& bm,
            const volVectorField& U
            // const volVectorField& Theta
        );

        //- Return the contribution to the linear momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> linearMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U
            // const volVectorField& Theta
        );

        //- Return the contribution to the angular momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> angularMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U
            // const volVectorField& Accl
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace momentumContributions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
