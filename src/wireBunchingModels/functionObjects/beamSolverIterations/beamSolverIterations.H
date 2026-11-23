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
    Foam::functionObjects::beamSolverIterations

Description
    Writes the solver work of a beam in each time step. One line per time
    step in postProcessing/<startTime>/<name>.dat:

        Time  nEvolve  nNewton  totalNewton

    nEvolve      calls to the beam's evolve() in the time step. In a coupled
                 run there is one per PIMPLE outer corrector that moves the
                 mesh; in a beamFoam run it is 1
    nNewton      Newton iterations (linear solves) in the time step, summed
                 over its evolve() calls
    totalNewton  nNewton summed since this run started

    In a time step in which the beam is not solved, nEvolve and nNewton
    are 0. Give the beam region when it is not the default region, e.g. in
    a coupled run:

    \verbatim
    beamSolverIterations
    {
        type    beamSolverIterations;
        libs    ("libwireBunchingModels.so");
        region  beamone;
    }
    \endverbatim

Author
    Colm Mc Alister, UCD.

SourceFiles
    beamSolverIterations.C

\*---------------------------------------------------------------------------*/

#ifndef beamSolverIterations_H
#define beamSolverIterations_H

#include "fvMeshFunctionObject.H"
#include "OFstream.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class beamModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class beamSolverIterations Declaration
\*---------------------------------------------------------------------------*/

class beamSolverIterations
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Reference to main object registry
        const Time& time_;

        //- Beam region name
        word regionName_;

        //- History file
        autoPtr<OFstream> historyFilePtr_;


    // Private Member Functions

        //- Write data
        bool writeData();


public:

    //- Runtime type information
    TypeName("beamSolverIterations");


    // Constructors

        //- Construct from Time and dictionary
        beamSolverIterations
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        beamSolverIterations(const beamSolverIterations&) = delete;

        //- No copy assignment
        void operator=(const beamSolverIterations&) = delete;


    //- Destructor
    virtual ~beamSolverIterations() = default;


    // Member Functions

        //- Read the region name
        virtual bool read(const dictionary& dict);

        //- Write the counts of the time step just solved
        virtual bool execute();

        //- Does nothing
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
