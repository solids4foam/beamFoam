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
    Foam::morisonForce

Description


SourceFiles
    morisonForce.C

Author
    Niels Gjoel Jacobsen, Deltares.  All rights reserved.


\*---------------------------------------------------------------------------*/

#ifndef morisonForce_H
#define morisonForce_H

#include "fvMesh.H"
#include "fvCFD.H"
#include "volFields.H"
#include "IOdictionary.H"

#include "autoPtr.H"
#include "runTimeSelectionTables.H"

//#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class morisonForce Declaration
\*---------------------------------------------------------------------------*/

class morisonForce
:
    public IOdictionary
{

protected:

    // Protected data
        const Time& runTime_;
        const fvMesh& mesh_;

        dimensionedScalar dragCoeff_;
        dimensionedScalar inertiaCoeff_;
        dimensionedScalar frictionCoeff_;

        dimensionedScalar rhoFluid_;
        dimensionedScalar rhoStem_;

        scalar G_;
        vector vertDir_;

        mutable volSymmTensorField addedMassCoeffs_;
        mutable volSymmTensorField dragCoeffs_;

        Switch linearize_;

    // Protected member functions
        inline scalar rampfactor(const scalar Tsoft) const
        {
            scalar res(1);

            if (SMALL < Tsoft)
            {
                scalar omega(2.0*M_PI/Tsoft);
                res = Foam::sin(omega/4.0*Foam::min(runTime_.time().value(), Tsoft));
            }

            return res;
        }


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        morisonForce(const morisonForce&);

        //- Disallow default bitwise assignment
        void operator=(const morisonForce&);

public:

    //- Runtime type information
    TypeName("morisonForce");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            morisonForce,
            dictionary,
            (
                const Time& runTime, const fvMesh& mesh, const dictionary& solidDict
            ),
            (runTime, mesh, solidDict)
        );

    // Constructors

        //- Construct from components
        morisonForce
        (
            const Time&,
            const fvMesh&,
            const dictionary&
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<morisonForce> New
        (
            const Time& runTime,
            const fvMesh& mesh,
            const dictionary& solidDict
        );

        virtual void printCoeffs() = 0;

    // Destructor

        virtual ~morisonForce();


    // Member Functions

        virtual dimensionedVector U
        (
            const point&,
            const scalar&
        ) const = 0;

        virtual dimensionedVector acc
        (
            const point&,
            const scalar&
        ) const = 0;

        virtual scalar Uabsolute() const = 0;


        void updateForce
        (
            volVectorField&,
            const volVectorField&,
            const volVectorField&,
            const tmp<vectorField>,
            const scalar&,
            const scalar&,
            const volScalarField&
        ) const;

        const symmTensorField& addedMassCoeff() const
        {
            return addedMassCoeffs_.internalField();
        };

        const symmTensorField& dragCoeff() const
        {
            return dragCoeffs_.internalField();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
