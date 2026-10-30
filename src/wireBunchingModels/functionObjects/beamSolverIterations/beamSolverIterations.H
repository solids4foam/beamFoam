/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::beamSolverIterations

Description
    Writes the solver work of a beam in each time step, to compare the
    monolithic and partitioned rigid body couplings. One line per time
    step in postProcessing/<startTime>/<name>.dat:

        Time  nPimple  nBeamSolves  nNewton  totalNewton

    nPimple      calls to beam.evolve() in the time step. In a coupled
                 interFoam run (moveMeshOuterCorrectors yes) there is one
                 per PIMPLE outer corrector; in a beamFoam run it is 1
    nBeamSolves  beam Newton loops (evolveBeam). Equal to nPimple except
                 with beamFoam's partitioned rigidBodyEnd, where each
                 coupling iteration is a Newton loop
    nNewton      calls to BlockEigenSolverOF::solve, i.e. Newton
                 iterations, summed over all loops of the time step
    totalNewton  nNewton summed from the start of the run

    Put it in the case's system/controlDict; in a coupled run give the
    beam region:

    \verbatim
    beamSolverIterations
    {
        type    beamSolverIterations;
        libs    ("libwireBunchingModels.so");
        region  beamone;
    }
    \endverbatim

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

        //- The beam model on regionName_
        const beamModel& beam() const;

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
