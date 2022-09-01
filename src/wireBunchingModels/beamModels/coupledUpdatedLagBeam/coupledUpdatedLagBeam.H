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
    coupledUpdatedLagBeam

Description

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    coupledUpdatedLagBeam.C

\*---------------------------------------------------------------------------*/

#ifndef coupledUpdatedLagBeam_H
#define coupledUpdatedLagBeam_H

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
                        Class coupledUpdatedLagBeam Declaration
\*---------------------------------------------------------------------------*/

class coupledUpdatedLagBeam
:
    public beamModel
{
    // Private data

        //- Increment of mean line displacement field
        mutable volVectorField DW_;

        //- Total mean line displacement field
        volVectorField W_;
  
        //- Increment of cross-section rotation field
        volVectorField DTheta_;

        //- Cross-section rotation field
        volVectorField theta_;

        //- Distributed moment
        volVectorField m_;
  
        //- Cross-section force
        surfaceVectorField Q_;
        surfaceVectorField dQ_;
  
        //- Cross-section moment
        surfaceVectorField M_;
        surfaceVectorField dM_;
  
        //- Cross-section rotation tensor
        mutable surfaceTensorField DLambda_;

        //- Cross-section rotation tensor
        surfaceTensorField Lambda_;
  
        //- Axial stretch ration
        surfaceScalarField stretchRatio_;
  
        //- Coefficient in the force equilibrium equation
        surfaceTensorField CQDW_;

        //- Coefficient in the force equilibrium equation
        surfaceTensorField CQDTheta_;
        
        //- Coefficient in the moment equilibrium equation
        surfaceTensorField CMDTheta_;
        
        //- Coefficient in the moment equilibrium equation
        surfaceTensorField CMDTheta2_;

        //- Mean line tangent vector
        surfaceVectorField dRuDs_;

        //- Strain vector
        surfaceVectorField epsilon_;
  
        //- Point mesh
        pointMesh pMesh_;

        //- Mean line point displacement field
        pointVectorField pointDW_;

        //- Block vector field for displacement (first entry) and
        // rotation (second entry)
        volVector6Field DWDTheta_;

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
        coupledUpdatedLagBeam(const coupledUpdatedLagBeam&);

        //- Disallow default bitwise assignment
        void operator=(const coupledUpdatedLagBeam&);

public:

    //- Runtime type information
    TypeName("coupledUpdatedLagBeam");

    // Constructors

        //- Construct from components
        coupledUpdatedLagBeam
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~coupledUpdatedLagBeam()
        {}

    // Member Functions

        // Access

            //- Each beamModel must indicate whether W or DW is the primary
            //  solution variable
            virtual volVectorField& solutionW() const
            {
                // This model solves for DW
                return DW_;
            }
  
            //- Each beamModel must indicate whether Lambda or DLambda
            //  is the primary solution variable
            virtual surfaceTensorField& solutionLambda() const
            {
                // This model solves for DLambda
                return DLambda_;
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
