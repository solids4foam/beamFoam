/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::rigidBodyEnd

Description
    Rigid body pinned to the end patch of a beam, owned by beamFoam.

    Holds the body mass, inertia and state, and advances it with the
    Newmark-beta method. The beam end is attached to the body at a pin, so
    the beam passes a force, and a torque through the arm from the centre of
    mass to the attachment point.

    Rotation is modelled when a moment of inertia is given. The orientation
    at the new time is exp(theta) & Q0, with theta the rotation vector of the
    time step, integrated with the same Newmark parameters as translation:
        J alpha + omega x (J omega) + cr omega = externalMoment + r x beamForce
    Without a moment of inertia the body is a point mass at the attachment.

    Motion constraints are given by the caller (setConstraints) as
    projections onto the free translation and rotation directions, in global
    axes. Accelerations, and in the monolithic coupling the Newton
    increments, are projected onto them, so constrained components keep
    their initial (zero) values.

    The body is read from the "rigidBodyEnd" sub-dictionary of the beam model
    coefficients:

        rigidBodyEnd
        {
            patch               right;   // default: endPatchName
            mass                10;
            momentOfInertia     (0 0 0); // principal, body axes;
                                         // zero: no rotation
            orientation         (1 0 0 0 1 0 0 0 1); // initial, body to
                                         // global axes
            centreOfMass        (2 0 0); // default: attachment point
            velocity            (0 0 0); // initial
            angularVelocity     (0 0 0); // initial, global axes
            externalForce       (0 0 0); // constant applied force [N]
            externalMoment      (0 0 0); // constant applied moment [N m]
            linearDamping       0;       // viscous force -c v [N s/m]
            angularDamping      0;       // viscous moment -c omega [N m s]
            newmarkBeta         0.25;
            newmarkGamma        0.5;
            nCouplingIterations 1;       // partitioned coupling iterations
            couplingTolerance   1e-8;    // change relative to
                                         // max(step, total) motion
            couplingRelaxation  1;       // fixed under-relaxation
            jacobianCheck       0;       // monolithic: finite-difference
                                         // check in the first N time steps
        }

    Gravity acts on the body through the beam model's g, at the centre of
    mass.

    Author
    Colm Mc Alister, UCD.

SourceFiles
    rigidBodyEnd.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyEnd_H
#define rigidBodyEnd_H

#include "Time.H"
#include "dictionary.H"
#include "diagTensor.H"
#include "point.H"
#include "OFstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class rigidBodyEnd Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyEnd
{
    // Private data

        //- Reference to time
        const Time& runTime_;

        //- Name of the beam patch the body is attached to
        word patchName_;

        //- Body mass
        scalar mass_;

        //- Principal moments of inertia, in the initial global axes
        diagTensor momentOfInertia_;

        //- True when the body rotates (non-zero moment of inertia)
        bool rotation_;

        //- Centre of mass, if given; otherwise the attachment point
        bool centreOfMassGiven_;
        point centreOfMass_;

        //- Arm from the centre of mass to the attachment point, in the
        //  initial configuration (global axes) and in body axes
        vector arm0_;
        vector armBody_;

        //- Constant applied force and moment
        vector externalForce_;
        vector externalMoment_;

        //- Linear viscous damping of translation and rotation
        scalar linearDamping_;
        scalar angularDamping_;

        //- Gravitational acceleration
        vector gravity_;

        //- Newmark parameters
        scalar beta_;
        scalar gamma_;

        //- Partitioned coupling controls
        label nCouplingIterations_;
        scalar couplingTolerance_;
        scalar couplingRelaxation_;

        //- Number of time steps with a finite-difference Jacobian check
        label nJacobianChecks_;

        //- Projections onto the free translation and rotation directions
        //  (global axes); I when unconstrained
        tensor tConstraints_;
        tensor rConstraints_;

        //- Displacement, velocity and acceleration at the new time
        vector x_;
        vector v_;
        vector a_;

        //- Displacement, velocity and acceleration at the old time
        vector x0_;
        vector v0_;
        vector a0_;

        //- Orientation, angular velocity and acceleration at the new time
        tensor Q_;
        vector omega_;
        vector alpha_;

        //- Orientation, angular velocity and acceleration at the old time
        tensor Q0_;
        vector omega0_;
        vector alpha0_;

        //- Rotation vector of the current time step
        vector theta_;

        //- Force from the beam on the body at the latest beam solution
        vector beamForce_;

        //- Torque about the centre of mass from the beam force
        vector beamTorque_;

        //- Beam end displacement W when the body displacement is zero
        vector attachmentOffset_;

        //- Time index of the current step
        label timeIndex_;

        //- Coupling iterations used in the current step
        label nIterationsUsed_;

        //- History file
        mutable autoPtr<OFstream> historyFilePtr_;


    // Private Member Functions

        //- Open the history file and write its header
        void makeHistoryFile();

        //- Angular Newmark predictor, acceleration and velocity for a step
        //  rotation vector
        vector angularPredictor() const;
        vector angularAcceleration(const vector& theta) const;
        vector angularVelocity(const vector& theta) const;

        //- Velocity at the new time for a displacement
        vector velocityAt(const vector& x) const;

        //- Moment of inertia in global axes at a step rotation
        tensor inertia(const vector& theta) const;


public:

    //- Runtime type information
    TypeName("rigidBodyEnd");


    // Constructors

        //- Construct from dictionary
        rigidBodyEnd
        (
            const Time& runTime,
            const dictionary& dict,
            const word& defaultPatchName,
            const vector& gravity
        );

        //- Disallow copy
        rigidBodyEnd(const rigidBodyEnd&) = delete;

        //- Disallow assignment
        void operator=(const rigidBodyEnd&) = delete;


    //- Destructor
    virtual ~rigidBodyEnd() = default;


    // Member Functions

        // Access

            const word& patchName() const
            {
                return patchName_;
            }

            scalar mass() const
            {
                return mass_;
            }

            const vector& displacement() const
            {
                return x_;
            }

            const vector& velocity() const
            {
                return v_;
            }

            const vector& acceleration() const
            {
                return a_;
            }

            const vector& beamForce() const
            {
                return beamForce_;
            }

            bool rotation() const
            {
                return rotation_;
            }

            //- Rotation vector of the accepted time step
            const vector& stepRotation() const
            {
                return theta_;
            }

            const vector& angularVelocity() const
            {
                return omega_;
            }

            const vector& angularAcceleration() const
            {
                return alpha_;
            }

            //- Orientation at the new time: body to global axes
            const tensor& Q() const
            {
                return Q_;
            }

            const vector& beamTorque() const
            {
                return beamTorque_;
            }

            //- Set the applied force and moment (global axes, moment about
            //  the centre of mass), e.g. from a fluid solver
            void setExternalLoad(const vector& force, const vector& moment)
            {
                externalForce_ = force;
                externalMoment_ = moment;
            }

            label nCouplingIterations() const
            {
                return nCouplingIterations_;
            }

            scalar couplingTolerance() const
            {
                return couplingTolerance_;
            }

            scalar couplingRelaxation() const
            {
                return couplingRelaxation_;
            }

            const vector& oldDisplacement() const
            {
                return x0_;
            }

            label nJacobianChecks() const
            {
                return nJacobianChecks_;
            }

            const vector& attachmentOffset() const
            {
                return attachmentOffset_;
            }

            const tensor& tConstraints() const
            {
                return tConstraints_;
            }

            const tensor& rConstraints() const
            {
                return rConstraints_;
            }

            //- Set the translation and rotation constraint projections
            //  (global axes) and project the current state onto them
            void setConstraints
            (
                const tensor& tConstraints,
                const tensor& rConstraints
            );

            void setAttachmentOffset(const vector& offset)
            {
                attachmentOffset_ = offset;
            }

            //- Set the initial attachment point; gives the arm from the
            //  centre of mass (which defaults to this point)
            void setAttachmentPoint(const point& attachmentPoint);


        // Evolution

            //- Total force on the body for a given force from the beam
            vector totalForce(const vector& beamForce) const;

            //- Store the old-time state when a new time step starts.
            //  Safe to call more than once per time step.
            void newTimeStep();

            //- Newmark displacement at the new time for a given beam force
            vector newmarkDisplacement(const vector& beamForce) const;

            //- Newmark predictor: displacement at the new time for zero
            //  new-time acceleration
            vector predictor() const;

            //- m/(beta*deltaT^2): inertia force per unit displacement
            //  beyond the predictor
            scalar inertiaCoefficient() const;

            //- Translational residual (force balance)
            //      m a + c v - externalForce - m g - beamForce
            vector translationalResidual
            (
                const vector& x,
                const vector& beamForce
            ) const;

            //- Derivative of the translational residual's inertia and
            //  damping terms with respect to x
            scalar translationalJacobian() const;

            //- Orientation at the new time for a step rotation vector
            tensor orientation(const vector& theta) const;

            //- Arm from the centre of mass to the attachment point
            vector arm(const vector& theta) const;

            //- Displacement of the attachment point
            vector attachmentDisplacement
            (
                const vector& x,
                const vector& theta
            ) const;

            //- Rotational residual (moment balance about the centre of mass)
            //      J alpha + omega x (J omega) + cr omega - externalMoment
            //    - arm x beamForce
            vector rotationalResidual
            (
                const vector& theta,
                const vector& beamForce
            ) const;

            //- Derivative of the rotational residual's inertia, gyroscopic
            //  and damping terms with respect to theta
            tensor rotationalInertiaJacobian(const vector& theta) const;

            //- Left Jacobian of the rotation: a change dtheta of the step
            //  rotation vector turns the body by the spatial rotation
            //  vector rotationJacobian(theta) & dtheta
            tensor rotationJacobian(const vector& theta) const;

            //- Newmark step rotation for a given beam force, with the
            //  inertia, arm and angular velocity from a previous estimate
            vector newmarkRotation
            (
                const vector& beamForce,
                const vector& thetaEstimate
            ) const;

            //- Accept a displacement, step rotation and the beam force they
            //  produced: accelerations and velocities follow from the
            //  Newmark relations
            void accept
            (
                const vector& displacement,
                const vector& theta,
                const vector& beamForce,
                const label nIterations
            );


        // Write

            //- Write one line of the history file
            void write() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
