 /*---------------------------------------------------------------------------* \
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
    coupledTotalLagNewtonRaphsonBeam

Description

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    coupledTotalLagNewtonRaphsonBeam.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTotalLagNewtonRaphsonBeam_H
#define coupledTotalLagNewtonRaphsonBeam_H

#include "beamModel.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvCFD.H"
//#include "multibeamFvBlockMatrix.H"
//#include "plasticityStressResultantReturn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

/*---------------------------------------------------------------------------*\
              Class coupledTotalLagNewtonRaphsonBeam Declaration
\*---------------------------------------------------------------------------*/

class coupledTotalLagNewtonRaphsonBeam
:
    public beamModel
{
    // Private data

        //- Mean line displacement field
        mutable volVectorField W_;

        //- Mean line displacement field increment
        mutable volVectorField WIncrement_;

        //- Mean line velocity field
        mutable volVectorField U_;

        //- Mean line acceleration field
 	      mutable volVectorField Accl_;

        //- Mean line displacement field increment
        volVectorField DW_;

        //- Cross-section rotation field
        mutable volVectorField Theta_;

        //- Cross-section rotation field
        mutable volVectorField ThetaIncrement_;

        //- Cross-section angular velocity field
        mutable volVectorField Omega_;
        //- Cross-section angular acceleration field
        mutable volVectorField dotOmega_;

        //- Cross-section rotation field increment
        volVectorField DTheta_;

        //- Switch to updated Lagrangian formulation
        bool updatedLagrangian_;

        //- Switch to Kirchhoff beam
        bool kirchhoffBeam_;

        //- Mean line displacement field
        mutable volVectorField totW_;

        //- Cross-section force
        surfaceVectorField Q_;

        //- Explicit contribution to cross-section force
        surfaceVectorField explicitQ_;

        //- Cross-section axial force
        surfaceScalarField Qa_;

        //- Cross-section moment
        surfaceVectorField M_;

        //- Cross-section moment (material)
        surfaceVectorField Mref_;

        //- Explicit contribution to cross-section moment
        surfaceVectorField explicitM_;

        //- Explicit contribution to cross-section moment
        surfaceVectorField explicitMQ_;

        //- Cross-section rotation tensor
        mutable surfaceTensorField Lambda_;

        //- Cell-centre rotation matrix
        volTensorField RM_;

        // //- Cell-centre torsion error
        // volScalarField torsionError_;

        // //- Cell-centre mean line tangent
        // volVectorField dRdS_;

        // //- Cross-section rotation tensor
        // surfaceTensorField T_;

        //- Cross-section reference point displacement field
        // surfaceVectorField refWf_;

		//- Mean line reference displacement at face-centres
        surfaceVectorField refWf_;

        //- Mean line reference displacement at cell-centres
        volVectorField refW_;

        //- Cross-section rotation tensor
        surfaceTensorField refLambda_;

        //- Cell-centre rotation matrix
        volTensorField refRM_;

        //- Mid-line tangent vector
        surfaceVectorField refTangent_;

        //- Axial stretch ration
        surfaceScalarField stretchRatio_;

        //- Coefficient in the force equilibrium equation
        surfaceTensorField CQW_;

        //- Coefficient in the force equilibrium equation
        surfaceTensorField CQTheta_;

        //- Coefficient in the force equilibrium equation
        surfaceTensorField CQDTheta_;

        //- Coefficient in the moment equilibrium equation
        surfaceTensorField CMTheta_;

        //- Coefficient in the moment equilibrium equation
        surfaceTensorField CMTheta2_;

        //- Coefficient in the force equilibrium equation
        surfaceTensorField CMQW_;

        //- Coefficient in the force equilibrium equation
        surfaceTensorField CMQTheta_;

        //- Mean line tangent vector
        surfaceVectorField dR0Ds_;

        //- Strain vector (translational)
        surfaceVectorField Gamma_;

        //- Kirchoff transformation tensor
        surfaceTensorField kirchhoffTransTensor_;

        //- Strain vector (rotational)
        surfaceVectorField K_;

        //- Point mesh
        pointMesh pMesh_;

        //- Mean line point displacement field
        pointVectorField pointW_;

        //- Block vector field for displacement (first entry) and
        // rotation (second entry)
        // volVector6Field WTheta_;

        //- Total Iterations count for all load steps
        scalar totalIter_;

        //- Number of CVs
       // scalar nCV_;

        // //- Young modulus
        // dimensionedScalar E_;

        // //- Shear modulus
        // dimensionedScalar G_;

        // //- Density
        // dimensionedScalar rho_;

        // //- Beam cross-section area
        // dimensionedScalar A_;

        // //- Beam cross-section constant
        // dimensionedScalar I_;

        // //- Beam cross-section constant
        // dimensionedScalar J_;

        // //- Beam cross-section constant
        // dimensionedScalar EI_;

        // //- Beam cross-section constant
        // dimensionedScalar GJ_;

        // //- Beam cross-section constant
        // dimensionedScalar EA_;

        // //- Beam cross-section constant
        // dimensionedScalar GA_;

        //- Beam cross-section constant
        mutable surfaceTensorField CQ_;
        // dimensionedDiagTensor CQ_;

        //- Beam cross-section constant
        mutable surfaceTensorField CDQ_;

        //- Beam cross-section constant
        surfaceTensorField CM_;
        // dimensionedDiagTensor CM_;

        //- Beam cross-section constant
        surfaceTensorField CDM_;

        //- Beam cross-section constant
        surfaceTensorField CDMDGamma_;

        //- Beam cross-section constant
        surfaceTensorField CDQDK_;
        //   //- Beam cross-section constant
         //   dimensionedTensor CI_;

 	      //- Beam momentum inertial term constant
 	      dimensionedScalar ARho_;

         //- Beam cross-section constant
         //   dimensionedTensor CI_;
 	      dimensionedTensor CIRho_;

        // Fields related to Newmark-beta time integration scheme

 	      //- Boolean to activate Newmark-beta time-integration scheme
 	      bool newmark_;

 	      //   //- Boolean to activate the tangent space version
 	      //   //- in time integration Newmark scheme
 	      //   bool tangentSpace_;

 	      //- Newmark beta coefficient
 	      scalar betaN_;

 	      //- Newmark gamma coefficient
 	      scalar gammaN_;
        //- Beam cross-section constant
        //dimensionedTensor CI_;

    // Fields related to plasticity

        //- Is it elasto-plastic model active
        // bool plasticity_;

        //- Plastic strain vector
        // surfaceVectorField GammaP_;

        //- Plastic strain vector
        // surfaceVectorField KP_;

        //- Run-time selectable method to calculate DGammaP and DKP
        //- eg radial return
        // autoPtr<plasticityStressResultantReturn>
        //     plasticityStressResultantReturnPtr_;

        //- Beam cell-centre curvatures
        volScalarField curvature_;

    // Parallel performance data

        //- Total time spend on contact update
        scalar totalContactTime_;

        //- Total time spend on solution of equilibrium equation
        scalar totalSolutionTime_;

        //- Processor index
        volScalarField proc_;

    // Private Member Functions

        //- Apply line contact force
        // void applyLineContact(multibeamFvBlockMatrix& eqn);

        // //- Apply point contact force
        // void applyPointContact(multibeamFvBlockMatrix& eqn);

        // //- Apply conical pulleys contact force
        // void applyConicalPulleysContact(multibeamFvBlockMatrix& eqn);

        // //- Apply conical pulleys contact force
        // void applyConicalPulleysContactNew(multibeamFvBlockMatrix& eqn);

        // //- Apply conical pulleys contact force
        // void applyToroidalPulleysContact(fvBlockMatrix<vector6>& eqn);

        //- Disallow default bitwise copy construct
        coupledTotalLagNewtonRaphsonBeam
        (
            const coupledTotalLagNewtonRaphsonBeam&
        );

        //- Disallow default bitwise assignment
        void operator=(const coupledTotalLagNewtonRaphsonBeam&);

public:

    //- Runtime type information
    TypeName("coupledTotalLagNewtonRaphsonBeam");

    // Constructors

        //- Construct from components
        coupledTotalLagNewtonRaphsonBeam
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~coupledTotalLagNewtonRaphsonBeam()
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

            //- Each beamModel must indicate whether W or DW is the primary
            //  solution variable
            virtual volVectorField& solutionDW() const
            {
                // This model solves for W
                return WIncrement_;
            }

            //- Each beamModel must indicate whether Theta or DTheta
            //  is the primary solution variable
            virtual volVectorField& solutionTheta() const
            {
                // This model solves for Theta
                return Theta_;
            }

            //- Each beamModel must indicate whether Lambda or DLambda
            //  is operating rotation matrix
            virtual surfaceTensorField& solutionLambda() const
            {
                // This model works with Lambda
                return Lambda_;
            }

            //- Each beamModel must indicate whether Lambda or DLambda
            //  is operating rotation matrix
            virtual const volTensorField& solutionRMC() const
            {
                // This model works with Lambda
                return RM_;
            }

            // //- Cell centre beam mean line tangent
            // virtual const volVectorField& solutionCellTangent() const
            // {
            //     // This model works with Lambda
            //     return dRdS_;
            // }

            //- Return current beam points (orderd from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentBeamPoints
            (
                const label bI = 0
            ) const;

            //- Return current beam points (orderd from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentGlobalBeamPoints
            (
                const label bI = 0
            ) const;

            //- Return current beam mid-line tangents (orderd from
            //  start to end patch) (face centres in 3d representation)
            virtual tmp<vectorField> currentBeamTangents
            (
                const label bI = 0
            ) const;

            //- Return current beam mid-line tangents
            //  (orderd from start to end patch)
            //  (face centres in 3d representation)
            virtual tmp<vectorField> currentGlobalBeamTangents
            (
                const label bI = 0
            ) const;

            //- Return beam points and tangents
            virtual void currentGlobalBeamPointsAndTangents
            (
                const label bI,
                vectorField& points,
                vectorField& tangents
            ) const;

            //- Return beam points and tangents
            virtual void currentGlobalBeamPointsAndTangents
            (
                const label bI,
                const labelListList& procSendPoints,
                // const labelListList& procReceivePoints,
                vectorField& points,
                vectorField& tangents
            ) const;

            //-
            virtual tmp<vectorField> currentDisplacementIncrement() const;

            //-
            virtual tmp<tensorField> currentRotationIncrement() const;

        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual scalar evolve();

            //- Update total fields
            virtual void updateTotalFields();

            //- Write fields
            virtual void writeFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
