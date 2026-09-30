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

    Holds the body mass, inertia and translational state, and advances the
    state with the Newmark-beta method. The attachment point is the body's
    centre of mass (point mass), so the body carries no torque; rotation is
    added in a later phase.

    The body is read from the "rigidBodyEnd" sub-dictionary of the beam model
    coefficients:

        rigidBodyEnd
        {
            patch               right;   // default: endPatchName
            mass                10;
            momentOfInertia     (0 0 0); // stored, not used yet
            externalForce       (0 0 0); // constant applied force [N]
            newmarkBeta         0.25;
            newmarkGamma        0.5;
            nCouplingIterations 1;       // partitioned coupling iterations
            couplingTolerance   1e-8;    // displacement change relative to
                                         // max(step, total) displacement
            couplingRelaxation  1;       // fixed under-relaxation
            jacobianCheck       0;       // monolithic: finite-difference
                                         // check in the first N time steps
        }

    Gravity acts on the body through the beam model's g.

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

        //- Principal moments of inertia (not used yet)
        diagTensor momentOfInertia_;

        //- Constant applied force
        vector externalForce_;

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

        //- Displacement, velocity and acceleration at the new time
        vector x_;
        vector v_;
        vector a_;

        //- Displacement, velocity and acceleration at the old time
        vector x0_;
        vector v0_;
        vector a0_;

        //- Force from the beam on the body at the latest beam solution
        vector beamForce_;

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

            void setAttachmentOffset(const vector& offset)
            {
                attachmentOffset_ = offset;
            }


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

            //- Accept a displacement and the beam force it produced:
            //  acceleration and velocity follow from the Newmark relations
            void accept
            (
                const vector& displacement,
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
