/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    linearElasticMisesPlastic

Description

SourceFiles
    linearElasticMisesPlastic.C

Author
    Zeljko Tukovic, FSB Zagreb

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMisesPlastic_H
#define linearElasticMisesPlastic_H

#include "plasticityStressResultantReturn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class beamModel;

/*---------------------------------------------------------------------------*\
                         Class linearElasticMisesPlastic Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMisesPlastic
:
    public plasticityStressResultantReturn
{
    // Private data

        // Reference to constitutive model
        beamModel& beamModel_;

        //- Equivalent plastic strain
        FieldField<Field, scalar> eqEP_;
    
        //- Plastic strain
        FieldField<Field, vector> EP_;
    
        //- Consistency parameter
        FieldField<Field, scalar> gamma_;
    
        //- Return direction to the yield surface
        FieldField<Field, vector> plasticN_;
    
        //- Initial yield stress
        dimensionedScalar sigmaY0_;
    
        //- Strain hardening modulus
        dimensionedScalar HP_;
    
        // // Yield stress
        // surfaceScalarField sigmaY_;

        // // Incremental change in sigmaY
        // surfaceScalarField DSigmaY_;

        // Increment of plastic strain
        surfaceVectorField DGammaP_;
        surfaceVectorField DKP_;

        //- Elasto-plastic moduli
        surfaceTensorField DQDGamma_;
        surfaceTensorField DMDK_;
        surfaceTensorField DMDGamma_;
        surfaceTensorField DQDK_;

        //- Number of yealing faces
        label nYieldingFaces_;

        // // Equivalent plastic strain increment
        // surfaceScalarField DEpsilonPEq_;

        // Active yielding flag
        // 1.0 for active yielding
        // 0.0 otherwise
        volScalarField activeYield_;

        //- List of internal faces for which
        //  cross-section fields will be writen
        labelList writeFaceIDs_;

        //- Plasticity on/off
        bool plasticity_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticMisesPlastic(const linearElasticMisesPlastic&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticMisesPlastic&);

        //- Yield function
        tmp<scalarField> yieldFunction
        (
            const vectorField& Strial,
            const scalarField& eP
        ) const;

        //- Yield function
        scalar yieldFunction
        (
            const vector& Strial,
            const scalar& eP
        ) const;

        //- Yield function derivative
        scalar yieldFunctionDerivative
        (
            const vector& S,
            const scalar& eqEP,
            const scalar& gamma
        ) const;

        //- Calc equivalent stress
        scalar vonMisesStress
        (
            const vector& S
        ) const;

        //- Calc equivalent stress
        tmp<scalarField> vonMisesStress
        (
            const vectorField& S
        ) const;
    
        //- Modified stiffness
        tmp<diagTensorField> barC
        (
            const scalarField& eqEP,
            const scalarField& gamma
        ) const;

        //- Modified stiffness
        diagTensor barC
        (
            const scalar eqEP,
            const scalar gamma
        ) const;

        // //- Calculate plastic strain increment
        // void calcPlasticStrainIncrement
        // (
        //     scalar& DCP,
        //     scalar& DEP,
        //     vector& DEpsilonP,
        //     tensor& DSDEpsilon,
        //     const scalar N,
        //     const vector& M,
        //     const scalar NP,
        //     const scalar CP,
        //     const vector dEpsilonP,
        //     const scalar dEP,
        //     const scalar Ftrial
        // );

        // //- Calculate plastic strain increment in case of pure plasticity
        // void calcPerfectPlasticStrainIncrement
        // (
        //     scalar& DCP,
        //     scalar& DEP,
        //     vector& DEpsilonP,
        //     tensor& DSDEpsilon,
        //     const scalar N,
        //     const vector& M,
        //     const scalar NP,
        //     const scalar CP,
        //     const vector dEpsilonP,
        //     const scalar dEP,
        //     const scalar Ftrial
        // );
    
public:

    //- Runtime type information
    TypeName("linearElasticMisesPlastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElasticMisesPlastic
        (
            const word& name,
            beamModel& beamModel
        );


    // Destructor

        virtual ~linearElasticMisesPlastic();


    // Member Functions

        //- Correct the rheological model
        virtual void correct();

        //- Return slave friction value fraction
        virtual const surfaceVectorField& DGammaP() const
        {
            return DGammaP_;
        }

        virtual const surfaceVectorField& DKP() const
        {
            return DKP_;
        }
    
        //- Elasto-plastic moduli
        virtual const surfaceTensorField& DQDGamma() const
        {
            return DQDGamma_;
        }
    
        virtual const surfaceTensorField& DMDK() const
        {
            return DMDK_;
        }

        virtual const surfaceTensorField& DMDGamma() const
        {
            return DMDGamma_;
        }
    
        virtual const surfaceTensorField& DQDK() const
        {
            return DQDK_;
        }
    
        //- Return ...
        virtual void updateYieldStress();

        //- Write plasticity informations
        virtual void writeInfo() const;
    
        //- Write law dictionary
        virtual void writeDict(Ostream& os) const {}
    
        //- Write plasticity informations
        virtual void writeFields() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
