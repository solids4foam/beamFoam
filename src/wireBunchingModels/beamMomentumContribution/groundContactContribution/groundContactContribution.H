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
    Foam::beamMomentumContribution::groundContactContribution

Description
    Derived class to evaluate the ground contact force contribution
    due to seabed interactions -
    which will be added to the beam solver governing eqns
    * For now, the seabed contact is only allowed in z-direction

SourceFiles
    groundContactContribution.C

Author
    Seevani Bali, UCD.
    Philip Cardiff, UCD.
    Amirhossein Taran, UCD.


\*---------------------------------------------------------------------------*/

#ifndef groundContactContribution_H
#define groundContactContribution_H

#include "beamMomentumContribution.H"
#include "HermiteSpline.H"
// #include "beamModel.H"
#include "coupledTotalLagNewtonRaphsonBeam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// namespace beamMomentumContributions
// {

/*---------------------------------------------------------------------------*\
                           Class groundContactContribution Declaration
\*---------------------------------------------------------------------------*/

class groundContactContribution
:
    public beamMomentumContribution
{

private:

    //- beam momentum contributions sub dictionary
    dictionary beamMomentumContribDict_;

    //-  damping coefficient for spring/damper ground model
    scalar gDamping_;

    //- stiffness coefficient for spring/damper ground model
    scalar gStiffness_;

    //- ground coordinate for ground contact model
    scalar groundZ_;

public:

    //- Runtime type information
    TypeName("groundContact");

    // Constructors

        //- Contruct from components
        groundContactContribution
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        ~groundContactContribution()
        {}


    // Member Functions

        //- Return the contribution to the block diagonal of the Jacobian matrix
        //  Each block should be a 6x6 matrix
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        // bm: Reference to beamModel object

        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        );

        //- Return the contribution to the linear momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> linearMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
         );

        //- Return the contribution to the angular momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> angularMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U,
            const volVectorField& Accl
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// } // End namespace beamMomentumContributions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
