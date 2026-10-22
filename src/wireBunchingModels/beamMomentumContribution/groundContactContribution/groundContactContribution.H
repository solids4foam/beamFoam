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
    The ground is a plane. A beam CV cell-centre below it is pushed out
    along the plane normal by a penalty spring (kNormal) and a damper
    (cNormal, acting on rebound only), and its sliding along the beam is
    resisted by friction (kTangential, muFriction).

    The plane is given by groundPoint and groundNormal, or by groundZ as a
    short form for the horizontal plane z = groundZ with normal (0 0 1):
    \verbatim
    groundContact
    {
        beamMomentumContributionType groundContact;

        groundContactCoeffs
        {
            groundPoint    (0 0 -0.5);  // a point on the ground
            groundNormal   (0 0 1);     // from the ground into the fluid
            // groundZ     -0.5;        // short form of the two entries above
            kNormal        1e3;
            cNormal        1.0;
            kTangential    0.0;
            muFriction     0.0;
        }
    }
    \endverbatim

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
#include "plane.H"

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

    //-  Normal stiffness coefficient for spring/damper ground model
    scalar kNormal_;

    //- Normal damping coefficient for spring/damper ground model
    scalar cNormal_;

    //- Tangential stiffness coefficient for spring/damper ground model
    scalar kTangential_;

    //- Friction coefficient for ground contact model
    scalar muFriction_;

    //- Ground plane, with the normal pointing from the ground into the fluid
    plane groundPlane_;

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


    // Static Member Functions

        //- Read the ground plane from groundPoint and groundNormal, or from
        //  groundZ as a short form for the plane z = groundZ with normal
        //  (0 0 1)
        static plane readGroundPlane(const dictionary& dict);


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
