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
    Foam::beamMomentumContribution

Description
    Abstract base class for adding external linear and angular momentum
    contributions to the beam solver.

    This class broadly follows the morisonForce class of Niels Gjoel Jacobsen,
    Deltares.

SourceFiles
    beamMomentumContribution.C

Author
    Philip Cardiff, UCD.
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef beamMomentumContribution_H
#define beamMomentumContribution_H

#include "fvMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class beamMomentumContribution Declaration
\*---------------------------------------------------------------------------*/

class beamMomentumContribution
// :
// public IOdictionary
{

private:

    // Private data
    const word name_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    beamMomentumContribution(const beamMomentumContribution&);

    //- Disallow default bitwise assignment
    void operator=(const beamMomentumContribution&);

public:

    //- Runtime type information
    TypeName("beamMomentumContribution");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            beamMomentumContribution,
            dictionary,
            (
                // const Time& runTime
             const word name,
             const dictionary& dict
            ),
            (name, dict)
        );

    // Constructors

        //- Construct from components
        beamMomentumContribution
        (
            const word& name,
            const dictionary& dict
        );


    // Selectors

        static autoPtr<beamMomentumContribution> New
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~beamMomentumContribution();


    // Member Functions

        //- Return the contribution to the block diagonal of the Jacobian matrix
        //  Each block should be a 6x6 matrix
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const beamModel& bm,
            const volVectorField& U
            // const volVectorField& Accl
        ) = 0;

        //- Return the contribution to the linear momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> linearMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U
            // const volVectorField& Accl
        ) = 0;

        //- Return the contribution to the angular momentum source
        //  U     : the current cell-centred velocity field
        //  Accl : the current cell-centred acceleration field
        virtual tmp<vectorField> angularMomentumSource
        (
            const beamModel& bm,
            const volVectorField& U
            // const volVectorField& Accl
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
