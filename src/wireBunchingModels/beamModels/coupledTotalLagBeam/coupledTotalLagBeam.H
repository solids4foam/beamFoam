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
    coupledTotalLagBeam

Description

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    coupledTotalLagBeam.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTotalLagBeam_H
#define coupledTotalLagBeam_H

#include "beamModel.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "fvBlockMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

/*---------------------------------------------------------------------------*\
                        Class coupledTotalLagBeam Declaration
\*---------------------------------------------------------------------------*/

class coupledTotalLagBeam
:
    public beamModel
{
    // Private data

        //- Mean line displacement field
        mutable volVectorField W_;

        //- Cross-section rotation field
        volVectorField Theta_;

        //- Cross-section rotation field
        volVectorField refTheta_;
  
        // //- Distributed moment
        // volVectorField m_;
  
        //- Cross-section force
        surfaceVectorField Q_;
  
        //- Explicit contribution to cross-section force
        surfaceVectorField explicitQ_;
  
        //- Cross-section moment
        surfaceVectorField M_;
  
        //- Explicit contribution to cross-section moment
        surfaceVectorField explicitM_;
  
        //- Cross-section rotation tensor
        mutable surfaceTensorField Lambda_;

        //- Cross-section rotation tensor
        surfaceTensorField refLambda_;
  
        //- Axial stretch ration
        surfaceScalarField stretchRatio_;
  
        //- Coefficient in the force equilibrium equation
        surfaceTensorField CQW_;

        //- Coefficient in the force equilibrium equation
        surfaceTensorField CQTheta_;
        
        //- Coefficient in the moment equilibrium equation
        surfaceTensorField CMTheta_;
        
        //- Mean line tangent vector
        surfaceVectorField dRuDs_;

        //- Mean line tangent vector
        surfaceVectorField i_;
  
        //- Strain vector
        surfaceVectorField epsilon_;
  
        //- Point mesh
        pointMesh pMesh_;

        //- Mean line point displacement field
        pointVectorField pointW_;

        //- Block vector field for displacement (first entry) and
        // rotation (second entry)
        volVector6Field WTheta_;

        //- Young modulus
        dimensionedScalar E_;
  
        //- Shear modulus
        dimensionedScalar G_;

        //- Beam cross-section area
        dimensionedScalar A_;

        //- Beam cross-section constant
        dimensionedScalar I_;

        //- Beam cross-section constant
        dimensionedScalar J_;
  
        //- Beam cross-section constant
        dimensionedScalar EI_;
  
        //- Beam cross-section constant
        dimensionedScalar GJ_;
  
        //- Beam cross-section constant
        dimensionedScalar EA_;
  
        //- Beam cross-section constant
        dimensionedScalar GA_;
  
        //- Beam cross-section constant
        dimensionedTensor CQ_;
  
        //- Beam cross-section constant
        dimensionedTensor CM_;
  
    // Private Member Functions

        //- Disallow default bitwise copy construct
        coupledTotalLagBeam(const coupledTotalLagBeam&);

        //- Disallow default bitwise assignment
        void operator=(const coupledTotalLagBeam&);

public:

    //- Runtime type information
    TypeName("coupledTotalLagBeam");

    // Constructors

        //- Construct from components
        coupledTotalLagBeam
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~coupledTotalLagBeam()
        {}

    // Member Functions

        // Access

            //- Each beamModel must indicate whether W or DW is the primary
            //  solution variable
            virtual volVectorField& solutionW() const
            {
                // This model solves for W
                return W_;
            }
  
            //- Each beamModel must indicate whether Lambda or DLambda
            //  is operating rotation matrix variable
            virtual surfaceTensorField& solutionLambda() const
            {
                // This model works with Lambda
                return Lambda_;
            }
  
            //- Return current beam points (orderd from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentBeamPoints(const label bI = 0) const;

            //-
            virtual tmp<vectorField> currentDisplacementIncrement() const;

            //-
            virtual tmp<tensorField> currentRotationIncrement() const;
  
        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual scalar evolve();

            //- Update total fields
            virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
