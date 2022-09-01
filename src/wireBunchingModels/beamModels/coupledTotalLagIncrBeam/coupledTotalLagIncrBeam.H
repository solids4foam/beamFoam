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
    coupledTotalLagIncrBeam

Description

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    coupledTotalLagIncrBeam.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTotalLagIncrBeam_H
#define coupledTotalLagIncrBeam_H

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
                        Class coupledTotalLagIncrBeam Declaration
\*---------------------------------------------------------------------------*/

class coupledTotalLagIncrBeam
:
    public beamModel
{
    // Private data

        //- Mean line displacement field
        mutable volVectorField DW_;

        //- Mean line displacement field
        volVectorField W_;
  
        //- Cross-section rotation field
        volVectorField DTheta_;

        //- Cross-section rotation field
        volVectorField refTheta_;
  
        //- Distributed moment
        volVectorField m_;
  
        //- Cross-section force
        surfaceVectorField Q_;
  
        //- Explicit contribution to cross-section force
        surfaceVectorField explicitQ_;
  
        //- Cross-section moment
        surfaceVectorField M_;
  
        //- Explicit contribution to cross-section moment
        surfaceVectorField explicitM_;
  
        //- Cross-section rotation tensor
        mutable surfaceTensorField DLambda_;
  
        //- Cross-section rotation tensor
        surfaceTensorField DT_;
  
        //- Cross-section rotation tensor
        surfaceTensorField Lambda_;

        //- Cross-section rotation tensor
        surfaceTensorField refLambda_;
  
        //- Cross-section rotation tensor
        surfaceTensorField refT_;
  
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

        //- Mean line tangent vector
        surfaceVectorField i_;
  
        //- Strain vector
        surfaceVectorField gamma_;
  
        //- Strain vector
        surfaceVectorField kappa_;
  
        //- Point mesh
        pointMesh pMesh_;

        //- Mean line point displacement field
        pointVectorField pointW_;

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
        coupledTotalLagIncrBeam(const coupledTotalLagIncrBeam&);

        //- Disallow default bitwise assignment
        void operator=(const coupledTotalLagIncrBeam&);

public:

    //- Runtime type information
    TypeName("coupledTotalLagIncrBeam");

    // Constructors

        //- Construct from components
        coupledTotalLagIncrBeam
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~coupledTotalLagIncrBeam()
        {}

    // Member Functions

        // Access

            //- Each solidModel must indicate whether W or DW is the primary
            //  solution variable
            virtual volVectorField& solutionW() const
            {
                // This model solves for DW
                return DW_;
            }
  
            //- Each solidModel must indicate whether Lambda or DLambda
            //  is operating rotation matrix variable
            virtual surfaceTensorField& solutionLambda() const
            {
                // This model solves for DW
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
