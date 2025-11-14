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
    coupledTotalLagNewtonRaphsonBeamEvolve.C
    coupledTotalLagNewtonRaphsonBeamMatrix.C
    coupledTotalLagNewtonRaphsonBeamBoundaryConditions.C

\*---------------------------------------------------------------------------*/

#ifndef coupledTotalLagNewtonRaphsonBeam_H
#define coupledTotalLagNewtonRaphsonBeam_H

#include "beamModel.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "momentumContribution.H"

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

    //- Fields related to time integration schemes

        //- ddtScheme name from fvSchemes dictionary
        //- Used to check whether steadyState or not
        const word ddtSchemeName_;

        //- d2dt2Scheme name from fvSchemes dictionary
        const word d2dt2SchemeName_;

    //- Fields related to Newmark-beta time integration scheme

        //- Boolean to activate Newmark-beta time-integration scheme
        // bool newmark_;

        //- Newmark beta coefficient
        scalar betaN_;

        //- Newmark gamma coefficient
        scalar gammaN_;

        //- Run-time selectable momentum contributions
        autoPtr<momentumContribution> momentumContribPtr_;

        //- Fields related to drag forces due to Morrison's equation
        //- Boolean to activate drag forces
        bool dragActive_;

        //- Normal drag coefficient
        scalar Cdn_;

        //- tangential drag coefficient
        scalar Cdt_;

        // Ground contact related parameters.

        //- Boolean to check wether ground contact is active or not
        bool groundContactActive_;

        //-  damping coefficient for spring/damper ground model.
        scalar gDamping_;

        //- stiffness coefficient for spring/damper ground model.
        scalar gStiffness_;

        //- ground coordinate for ground contact Model.
        scalar groundZ_;

    // Parallel performance data

        //- Total time spend on contact update
        scalar totalContactTime_;

        //- Total time spend on solution of equilibrium equation
        scalar totalSolutionTime_;

        //- Processor index
        volScalarField proc_;

    // Private Member Functions

        //- Update the coefficients associated with the primary
        //- variables W_ and Theta_
        //- Function call implemented in
        //- 'coupledTotalLagNewtonRaphsonBeamEvolve.C'
        void updateEqnCoefficients();

        //- Assemble the l,d,u coefficients to solve for
        //- variables W_ and Theta_ in block coupled manner
        //- Function call found in 'coupledTotalLagNewtonRaphsonBeamMatrix.C'
        void assembleMatrixCoefficients
        (
            Field<scalarSquareMatrix>& d,
            Field<scalarSquareMatrix>& l,
            Field<scalarSquareMatrix>& u,
            Field<scalarRectangularMatrix>& source
        );

        //- Add the boundary contributions to l,d,u fields
        //- Member function found in coupledTotalLagNewtonRaphsonBeamMatrix.C
        //- This function for BC implementation can be found in
        //- 'coupledTotalLagNewtonRaphsonBeamBoundaryConditions.C'
        void assembleBoundaryConditions
        (
            Field<scalarSquareMatrix>& d,
            Field<scalarSquareMatrix>& l,
            Field<scalarSquareMatrix>& u,
            Field<scalarRectangularMatrix>& source
        );

        //- Update the solution and output variables like displacements,
        //- rotations, strains, forces and moments
        //- Function call implemented in
        //- 'coupledTotalLagNewtonRaphsonBeamEvolve.C'
        void updateSolutionVariables();

        //- Convergence check given all settings
        //  Determines whether the nonlinear solver has converged based on
        //  various criteria, including residual norms, step norms, and
        //  iteration limits.
        //  Similar approach to PETSc SNES solver.
        //  Arguments:
        //  - currentResidualNorm: The 2-norm of the current residual,
        //    indicating the deviation of the current solution from equilibrium.
        //  - initialResidualNorm: The 2-norm of the initial residual, used to
        //    evaluate relative convergence.
        //  - deltaXNorm: The 2-norm of the last step (increment in the
        //    solution vector) taken by the solver.
        //  - xNorm: The 2-norm of the current solution vector, used to scale
        //    the step tolerance for relative checks.
        //  - iteration: The current iteration count of the nonlinear solver.
        //  - maxIterations: The maximum number of allowed iterations before
        //    the solver terminates.
        //  - rtol: The relative residual norm tolerance. Convergence is
        //    achieved if the current residual is less than
        //    rtol*initialResidualNorm.
        //  - atol: The absolute residual norm tolerance. Convergence is
        //    achieved if the current residual is less than atol.
        //  - stol: The relative step norm tolerance. Convergence is achieved
        //    if deltaXNorm is less than stol*xNorm.
        //  - divtol: The divergence tolerance. If the current residual norm
        //    grows beyond divtol*initialResidualNorm, the solver is
        //    considered diverged.
        //  - writeResidualFrequency: Frequency at which the residuals are
        //    printed to the console. Set to 0 or a negative number to disable
        //    output.
        bool checkConvergence
        (
            const scalar currentResidualNorm,
            const scalar initialResidualNorm,
            const scalar deltaXNorm,
            const scalar xNorm,
            const label iteration,
            const label maxIterations,
            const scalar rtol,
            const scalar atol,
            const scalar stol,
            const scalar divtol,
            const label writeResidualFrequency,
            const bool writeConvergedReason = true
        );

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

        // Edit

            //- Evolve the beam solver and solve the mathematical model
            // - Function call in 'coupledTotalLagNewtonRaphsonBeamEvolve.C'
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
