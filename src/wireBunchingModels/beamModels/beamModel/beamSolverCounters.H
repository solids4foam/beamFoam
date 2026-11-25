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
    Foam::beamSolverCounters

Description
    Counts the work of a beam solver: calls to evolve(), beam solves (beam
    Newton loops) and Newton iterations (linear solves) in the current time
    step, and Newton iterations since this run started (a restarted run
    counts from zero). A beam model updates them, and the
    beamConvergenceData function object writes them. The count of beam
    solves in the time step also tells the beam model when a time step
    starts (e.g. for the Newmark predictor), so it is part of the solver
    state within a time step.

    The counts of a time step start at the first evolve() of that time
    step. A coupled solver calls evolve() once per PIMPLE outer corrector
    that moves the mesh; beamFoam calls it once per time step. Each
    evolve() solves the beam once, or once per coupling iteration with a
    partitioned rigidBodyEnd.

Author
    Colm Mc Alister, UCD.

\*---------------------------------------------------------------------------*/

#ifndef beamSolverCounters_H
#define beamSolverCounters_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class beamSolverCounters Declaration
\*---------------------------------------------------------------------------*/

class beamSolverCounters
{
    // Private Data

        //- Time index of the per time step counts
        label timeIndex_;

        //- Calls to evolve() in time step timeIndex_
        label nEvolveCalls_;

        //- Beam solves (beam Newton loops) in time step timeIndex_
        label nBeamSolves_;

        //- Newton iterations in time step timeIndex_
        label nNewtonIterations_;

        //- Newton iterations since this run started
        label totalNewtonIterations_;


public:

    // Constructors

        //- Construct with all counts zero
        beamSolverCounters()
        :
            timeIndex_(-1),
            nEvolveCalls_(0),
            nBeamSolves_(0),
            nNewtonIterations_(0),
            totalNewtonIterations_(0)
        {}


    // Member Functions

        //- Count a call to evolve() in time step timeIndex. The first call
        //  of a time step resets the per time step counts
        void countEvolveCall(const label timeIndex)
        {
            if (timeIndex != timeIndex_)
            {
                timeIndex_ = timeIndex;
                nEvolveCalls_ = 0;
                nBeamSolves_ = 0;
                nNewtonIterations_ = 0;
            }

            ++nEvolveCalls_;
        }

        //- Count a beam solve (beam Newton loop) in the time step of the
        //  latest countEvolveCall()
        void countBeamSolve()
        {
            ++nBeamSolves_;
        }

        //- Add the Newton iterations (linear solves) of one beam solve,
        //  i.e. the beam's iOuterCorr() at the end of the solve
        void addNewtonIterations(const label nIter)
        {
            nNewtonIterations_ += nIter;
            totalNewtonIterations_ += nIter;
        }

        //- Calls to evolve() in time step timeIndex: zero if the beam was
        //  not solved in that time step
        label nEvolveCalls(const label timeIndex) const
        {
            return timeIndex == timeIndex_ ? nEvolveCalls_ : 0;
        }

        //- Beam solves in time step timeIndex: zero if the beam was not
        //  solved in that time step
        label nBeamSolves(const label timeIndex) const
        {
            return timeIndex == timeIndex_ ? nBeamSolves_ : 0;
        }

        //- Newton iterations in time step timeIndex: zero if the beam was
        //  not solved in that time step
        label nNewtonIterations(const label timeIndex) const
        {
            return timeIndex == timeIndex_ ? nNewtonIterations_ : 0;
        }

        //- Newton iterations since this run started
        label totalNewtonIterations() const
        {
            return totalNewtonIterations_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
