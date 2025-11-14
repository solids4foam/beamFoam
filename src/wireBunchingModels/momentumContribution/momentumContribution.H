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
    Foam::momentumContribution

Description
    Abstract base class for adding external linear and angular momentum
    contributions to the beam solver.

    This class broadly follows the morisonForce class of Niels Gjoel Jacobsen,
    Deltares.

SourceFiles
    momentumContribution.C

Author
    Philip Cardiff, UCD.
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef momentumContribution_H
#define momentumContribution_H

#include "fvMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class momentumContribution Declaration
\*---------------------------------------------------------------------------*/

class momentumContribution
:
    public IOdictionary
{

private:

    // Private data

        //- Read-access to the mesh
        const Time& runTime_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        momentumContribution(const momentumContribution&);

        //- Disallow default bitwise assignment
        void operator=(const momentumContribution&);

public:

    //- Runtime type information
    TypeName("momentumContribution");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            momentumContribution,
            dictionary,
            (
                const Time& runTime
            ),
            (runTime)
        );

    // Constructors

        //- Construct from components
        momentumContribution
        (
            const Time&
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<momentumContribution> New
        (
            const Time& runTime
        );


    // Destructor

        virtual ~momentumContribution();


    // Member Functions

        //- Return the contribution to the block diagonal of the Jacobian matrix
        //  Each block should be a 6x6 matrix
        //  W     : the current cell-centred displacement field
        //  Theta : the current cell-centred rotation field
        virtual tmp<Field<scalarSquareMatrix>> diagCoeff
        (
            const volVectorField& W, const volVectorField& Theta
        ) = 0;

        //- Return the contribution to the linear momentum source
        //  W     : the current cell-centred displacement field
        //  Theta : the current cell-centred rotation field
        virtual tmp<vectorField> linearMomentumSource
        (
            const volVectorField& W, const volVectorField& Theta
        ) = 0;

        //- Return the contribution to the angular momentum source
        //  W     : the current cell-centred displacement field
        //  Theta : the current cell-centred rotation field
        virtual tmp<vectorField> angularMomentumSource
        (
            const volVectorField& W, const volVectorField& Theta
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
