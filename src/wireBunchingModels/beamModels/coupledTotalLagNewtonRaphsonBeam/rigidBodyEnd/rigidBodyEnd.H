/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::rigidBodyEnd

Description
    Rigid body attached to the end patch of a beam.

    Holds the body mass, inertia and state, and advances it with the
    Newmark-beta method. The joint between the beam end and the body is
    either
        pin      the beam passes a force, and a torque through the arm from
                 the centre of mass to the attachment point (default). The
                 beam end must be free to rotate
        welded   the beam end cross-section also turns with the body, and
                 the beam passes its end moment as well; needs a moment of
                 inertia and a fixedValue Theta condition on the patch

    Rotation is modelled when a moment of inertia is given. The orientation
    at the new time is rotationMatrix(theta) & Q0, with theta the rotation
    vector of the time step, integrated with the same Newmark parameters as
    translation:
        J alpha + omega x (J omega) + cr omega
          = externalMoment + r x beamForce + beamMoment
    with beamMoment = 0 for a pin.
    Without a moment of inertia the body is a point mass at the attachment.

    The body and the beam are solved in turn (partitioned coupling): each
    coupling iteration moves the body with the beam force and moment of the
    previous beam solve, and solves the beam with the end set by the body.
    A time step stops iterating when the body motion changes by less than
    couplingTolerance, relative to the larger of its motion in the step and
    its total motion, or after nCouplingIterations beam solves. A time step
    that reaches the limit gives a warning and is counted in the history
    file.

    Stability: with a single coupling iteration (staggered coupling) the body
    only sees the beam force of the previous solve. For a welded body twisted
    by the beam, beta deltaT^2 kt/J (kt the torsional stiffness of the beam
    at the end, J the body's moment of inertia) indicates how hard the
    coupling is: well below 1 the iterations converge quickly, and around 1
    and above they slow down and may diverge. This is a guide, not a bound:
    other modes (axial, bending, a pin with an offset centre of mass) have
    their own stiffness-to-inertia ratios. Reduce deltaT or use under-
    relaxation when the coupling does not converge.

    The body state is written to <time>/uniform/rigidBodyEnd at every write
    time (of the beam region) and read back when the run restarts from that
    time. One line per time step is written to
    postProcessing/rigidBodyEnd/<startTime>/rigidBodyEnd.dat, or
    postProcessing/<region>/rigidBodyEnd/<startTime>/rigidBodyEnd.dat for a
    beam that is not in the default region (as for function objects).

    The body is read from the "rigidBodyEnd" sub-dictionary of the beam model
    coefficients:

        rigidBodyEnd
        {
            patch               right;   // default: endPatchName
            joint               pin;     // pin or welded
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
            nCouplingIterations 20;      // beam solves per time step
            couplingTolerance   1e-8;    // change relative to
                                         // max(step, total) motion
            couplingRelaxation  1;       // fixed under-relaxation
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
#include "IOdictionary.H"
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
    // Private Data

        //- Reference to time
        const Time& runTime_;

        //- Name of the beam patch the body is attached to
        word patchName_;

        //- Joint between the beam end and the body
        enum class jointType { pin, welded };
        jointType joint_;

        //- Body mass
        scalar mass_;

        //- Principal moments of inertia, in body axes
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

        //- Body state, written at write times and read on restart
        IOdictionary state_;

        //- True when the state was read from a previous run
        bool restarted_;

        //- Displacement, velocity and acceleration at the new time
        vector x_;
        vector v_;
        vector a_;

        //- Displacement, velocity and acceleration at the old time
        vector x0_;
        vector v0_;
        vector a0_;

        //- Orientation (body to global axes), angular velocity and
        //  acceleration at the new time
        tensor Q_;
        vector omega_;
        vector alpha_;

        //- Initial orientation
        tensor Qinit_;

        //- Orientation, angular velocity and acceleration at the old time
        tensor Q0_;
        vector omega0_;
        vector alpha0_;

        //- Rotation vector of the current time step
        vector theta_;

        //- Sum of the step rotation vectors since the start, at the new and
        //  the old time. For rotation about a fixed axis it is the rotation
        //  angle times the axis, and unlike a rotation vector of Q it does
        //  not wrap at pi
        vector rotationSum_;
        vector rotationSum0_;

        //- Force from the beam on the body at the latest beam solution
        vector beamForce_;

        //- Moment from the beam on the body (welded joint; zero for a pin)
        vector beamMoment_;

        //- Total torque from the beam about the centre of mass:
        //  arm x beamForce + beamMoment
        vector beamTorque_;

        //- Beam end displacement W when the body displacement is zero
        vector attachmentOffset_;

        //- Time index of the current step
        label timeIndex_;

        //- Coupling iterations (beam solves) used in the current step
        label nIterationsUsed_;

        //- True when the current step's coupling iterations converged
        bool converged_;

        //- Time steps whose coupling iterations did not converge, since the
        //  start of the run (including before a restart)
        label nNotConverged_;

        //- History file
        mutable autoPtr<OFstream> historyFilePtr_;


    // Private Member Functions

        //- Read the state written by a previous run
        void readState();

        //- Copy the state into state_, for writing
        void storeState();

        //- Open the history file of the body on the mesh region db and
        //  write its header
        void makeHistoryFile(const objectRegistry& db);

        //- Newmark predictor: displacement at the new time for zero
        //  new-time acceleration
        vector predictor() const;

        //- Angular Newmark predictor, acceleration and velocity for a step
        //  rotation vector
        vector angularPredictor() const;
        vector angularAcceleration(const vector& theta) const;
        vector angularVelocity(const vector& theta) const;

        //- Moment of inertia in global axes at a step rotation
        tensor inertia(const vector& theta) const;

        //- Orientation at the new time for a step rotation vector
        tensor orientation(const vector& theta) const;

        //- Total force on the body for a given force from the beam
        vector totalForce(const vector& beamForce) const;


public:

    //- Runtime type information
    TypeName("rigidBodyEnd");


    // Constructors

        //- Construct from dictionary. The state is read from
        //  <time>/uniform/rigidBodyEnd of the registry (the beam mesh) when
        //  it exists
        rigidBodyEnd
        (
            const objectRegistry& db,
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

            //- True when the body rotates
            bool rotation() const
            {
                return rotation_;
            }

            //- True for a welded joint: the beam end turns with the body
            bool welded() const
            {
                return joint_ == jointType::welded;
            }

            //- True when the state was read from a previous run
            bool restarted() const
            {
                return restarted_;
            }

            //- Radius of gyration about the largest principal axis
            scalar radiusOfGyration() const
            {
                return Foam::sqrt(cmptMax(momentOfInertia_)/mass_);
            }

            const vector& displacement() const
            {
                return x_;
            }

            const vector& oldDisplacement() const
            {
                return x0_;
            }

            const vector& velocity() const
            {
                return v_;
            }

            const vector& angularVelocity() const
            {
                return omega_;
            }

            const vector& beamForce() const
            {
                return beamForce_;
            }

            const vector& beamMoment() const
            {
                return beamMoment_;
            }

            const vector& attachmentOffset() const
            {
                return attachmentOffset_;
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


        // Set-up

            //- Attach the body to the beam end. endPoint is the end face
            //  centre without the beam displacement W, and endW the current
            //  W there. At the start of a run the current position is the
            //  attachment point, and the accelerations follow from the loads
            //  including the current force and moment from the beam on the
            //  body; on a restart both come from the saved state
            void attach
            (
                const point& endPoint,
                const vector& endW,
                const vector& beamForce,
                const vector& beamMoment
            );


        // Evolution

            //- Store the old-time state when a new time step starts.
            //  Safe to call more than once per time step.
            void newTimeStep();

            //- Newmark displacement at the new time for a given beam force
            vector newmarkDisplacement(const vector& beamForce) const;

            //- Newmark step rotation for a given beam force and moment, with
            //  the inertia, arm and angular velocity from a previous estimate
            vector newmarkRotation
            (
                const vector& beamForce,
                const vector& thetaEstimate,
                const vector& beamMoment
            ) const;

            //- Angle turned since the start, for a step rotation vector: the
            //  magnitude of the sum of the step rotation vectors. Does not
            //  wrap at pi
            scalar totalRotation(const vector& theta) const;

            //- Arm from the centre of mass to the attachment point
            vector arm(const vector& theta) const;

            //- Displacement of the attachment point
            vector attachmentDisplacement
            (
                const vector& x,
                const vector& theta
            ) const;

            //- Rotation of the body over the time step, for a step rotation
            //  vector
            tensor stepRotation(const vector& theta) const;

            //- Accept a displacement, step rotation and the beam force and
            //  moment they produced: accelerations and velocities follow
            //  from the Newmark relations
            void accept
            (
                const vector& displacement,
                const vector& theta,
                const vector& beamForce,
                const vector& beamMoment,
                const label nIterations,
                const bool converged
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
