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
    Foam::functionObjects::beamConvergenceData

Description
    Writes the convergence and solver work of a beam in each time step.
    One line per time step in postProcessing/<startTime>/
    beamConvergenceData.dat:

        Time  nCorrectors  nEvolve  nNewton  nBeamSolves

    nCorrectors  Newton iterations in the last beam solve (beam Newton
                 loop) in the time step
    nEvolve      calls to the beam's evolve() in the time step. In a coupled
                 run there is one per PIMPLE outer corrector that moves the
                 mesh; in a beamFoam run it is 1
    nNewton      Newton iterations (linear solves) in the time step, summed
                 over its beam solves. Equal to nCorrectors when
                 nBeamSolves is 1
    nBeamSolves  beam Newton loops in the time step: one per evolve(), or
                 one per coupling iteration with a partitioned rigidBodyEnd

    In a time step in which the beam is not solved, nEvolve, nNewton and
    nBeamSolves are 0. The total Newton iterations of a run is the sum of nNewton. Give
    the beam region when it is not the default region, e.g. in a coupled
    run:

    \verbatim
    beamConvergenceData
    {
        type    beamConvergenceData;
        libs    ("libwireBunchingModels.so");
        region  beamone;
    }
    \endverbatim

Author
    Seevani Bali, UCD.  All rights reserved
    Colm Mc Alister, UCD. Solver work columns (nEvolve, nNewton,
    nBeamSolves)

SourceFiles
    beamConvergenceData.C

\*---------------------------------------------------------------------------*/

#ifndef beamConvergenceData_H
#define beamConvergenceData_H

#include "fvMeshFunctionObject.H"
#include "OFstream.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class beamConvergenceData Declaration
\*---------------------------------------------------------------------------*/

class beamConvergenceData
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Forces/moment file ptr
        autoPtr<OFstream> historyFilePtr_;

    // Private Member Functions

        //- Write data
        bool writeData();

public:

    //- Runtime type information
    TypeName("beamConvergenceData");


    // Constructors

        //- Construct from Time and dictionary
        beamConvergenceData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        beamConvergenceData(const beamConvergenceData&) = delete;

        //- No copy assignment
        void operator=(const beamConvergenceData&) = delete;


    //- Destructor
    virtual ~beamConvergenceData() = default;


    // Member Functions

        //- Read the beamConvergenceData data
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the beamConvergenceData
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
