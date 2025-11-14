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
    Foam::momentumContributions::morisonMomentum

Description
    SEEVANI!!!

SourceFiles
    morisonMomentum.C

Author
    Seevani Balii, UCD.
    Philip Cradiff, UCD.


\*---------------------------------------------------------------------------*/

#ifndef morisonMomentum_H
#define morisonMomentum_H

#include "momentumContribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace momentumContributions
{

/*---------------------------------------------------------------------------*\
                           Class morisonMomentum Declaration
\*---------------------------------------------------------------------------*/

class morisonMomentum
:
    public momentumContribution
{

private:

    // Protected data

        //- SEEVANI
        // vector Um_;
        // scalar Tsoft_;

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
        //  W     : the current cell-centred displacement field
        //  Theta : the current cell-centred rotation field
        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const volVectorField& W, const volVectorField& Theta
        );

        //- Return the contribution to the linear momentum source
        //  W     : the current cell-centred displacement field
        //  Theta : the current cell-centred rotation field
        virtual tmp<vectorField> linearMomentumSource
        (
            const volVectorField& W, const volVectorField& Theta
        );

        //- Return the contribution to the angular momentum source
        //  W     : the current cell-centred displacement field
        //  Theta : the current cell-centred rotation field
        virtual tmp<vectorField> angularMomentumSource
        (
            const volVectorField& W, const volVectorField& Theta
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace momentumContributions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
