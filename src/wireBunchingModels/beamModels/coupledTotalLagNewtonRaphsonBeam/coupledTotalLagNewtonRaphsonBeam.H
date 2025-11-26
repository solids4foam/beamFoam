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
    This code is the total Lagrangian formulation to solve displacements (W)
    and rotations (Theta) implemented in a coupled Newton-Raphson approach for
    geometrically exact Simo-Reissner beams

    More details on the implementation can be found in the paper:
    https://doi.org/10.1002/nme.6994

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Seevani Bali, UCD.

SourceFiles
    coupledTotalLagNewtonRaphsonBeam.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTotalLagNewtonRaphsonBeam_H
#define coupledTotalLagNewtonRaphsonBeam_H

#include "beamModel.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "meshSearch.H"
#include "vectorIOList.H"
#include "vectorList.H"
//#include "multibeamFvBlockMatrix.H"

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

        //- Fluid Velocity around the beam
        mutable volVectorField fluidU_;

        //- Fluid cell index that has been found for each beam cell
        mutable labelIOList fluidCellIDs_;

        //- immersed force acting on each beam cell
        mutable vectorIOList immersedForce_;

        //- Cell marker function for drag calculation
        mutable volScalarField cellMarker_;

        //- beam cell center veclocity stored on fluid mesh
        mutable volVectorField beamVelocity_;

        //- actuator line method force acting on each beam cell
        mutable volVectorField almForce_;

        //- Mean line acceleration field
 	    mutable volVectorField Accl_;

        //- Mean line displacement field increment
        volVectorField DW_;

        //- Cross-section rotation field
        mutable volVectorField Theta_;

        //- Cross-section angular velocity field
        mutable volVectorField Omega_;

        //- Cross-section angular acceleration field
        mutable volVectorField dotOmega_;

        //- Cross-section rotation field increment
        volVectorField DTheta_;

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

        //- Face-centre rotation tensor field
        mutable surfaceTensorField Lambdaf_;

        //- Cell-centre rotation tensor field
        volTensorField Lambda_;

        // //- Cell-centre mean line tangent
        // volVectorField dRdS_;

		//- Mean line reference displacement at face-centres
        surfaceVectorField refWf_;

        //- Mean line reference displacement at cell-centres
        volVectorField refW_;

        //- Face-centre reference rotation tensor field
        surfaceTensorField refLambdaf_;

        //- Cell-centre reference rotation tensor field
        volTensorField refLambda_;

        //- Mid-line tangent vector field at face centres
        surfaceVectorField refTangent_;

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

        //- Beam cross-section constant
        mutable surfaceTensorField CQ_;

        //- Beam cross-section constant
        surfaceTensorField CM_;

        //- Beam cross-section constant
        surfaceTensorField CDQDK_;

 	    //- Beam momentum inertial term constant
 	    dimensionedScalar ARho_;

        //- Beam cross-section constant
 	    dimensionedTensor CIRho_;

        // Fields related to Newmark-beta time integration scheme

        //- Boolean to activate Newmark-beta time-integration scheme
        bool newmark_;

        //- Newmark beta coefficient
        scalar betaN_;

        //- Newmark gamma coefficient
        scalar gammaN_;
        //- Fields related to drag forces due to Morrison's equation
        //- Boolean to activate drag forces
        bool dragActive_;

        //- Boolean to enable beam-fluid coupling and its effect on drag calculations
        bool beamFluidInteraction_;

        //- Boolean to enable immersed force on beam cells
        bool immersedForceActive_;

        //- Boolean to enable ALM like sampling of the fluid velocity per each beam cell
        bool almSamplingActive_;

        //- Boolean to enable Morison Force drag and added mass force
        bool MorisonForceActive_;

        //- Normal drag coefficient
        scalar Cdn_;

        //- tangential drag coefficient
        scalar Cdt_;

        //- Normal added mass coefficient
        scalar CMn_;

        //- Tangential added mass coefficient
        scalar CMt_;

        //- Seabed friction parameters

        //- Boolean to check wether seabed friction is active or not
        bool seadbedFrictionActive_;

        //- seabed friction coefficient
        scalar muFriction_;


        // Ground contact related parameters.

        //- Boolean to check wether ground contact is active or not
        bool groundContactActive_;

        //-  damping coefficient for spring/damper ground model
        scalar cNormal_;

        //- stiffness coefficient for spring/damper ground model, in normal direction
        scalar kNormal_;

        //- stiffness coefficient for seabed friction model, in tangential direction
        scalar kTangential_;

        //- ground coordinate for ground contact Model.
        scalar groundZ_;

        //- meshSearch engine for finding adjacent fluid Cells
        meshSearch searchEngine_;

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
                // This model works with Lambdaf
                return Lambdaf_;
            }

            //- Each beamModel must indicate whether Lambda or DLambda
            //  is operating rotation matrix
            virtual const volTensorField& solutionRMC() const
            {
                // This model works with Lambda
                return Lambda_;
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
;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
