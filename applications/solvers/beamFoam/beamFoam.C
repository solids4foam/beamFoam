/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

Application
    beamFoam

Description
    Solves a Timoshenko beam models for
    non-linear elastic beams (small strains and large rotations).

\*---------------------------------------------------------------------------*/

// #include "mpi.h"
#include "parRun.H"

#include "objectRegistry.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "linear.H"
#include "uniformDimensionedFields.H"
#include "calculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "mathematicalConstants.H"

#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"

// #include <Eigen/Core>

// #ifndef namespaceFoam
// #define namespaceFoam
//     using namespace Foam;
// #endif


#include "fvCFD.H"
#include "beamModel.H"

// #include "petsc.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Eigen::initParallel();
    // Eigen::setNbThreads(2);
    // Foam::Info << "nThreads: " << Eigen::nbThreads() << Foam::endl;

    // PetscErrorCode ierr =
    //     PetscInitializeNoArguments();
    // if (ierr)
    // {
    //     return ierr;
    // }

    // MPI_Init(&argc, &argv);

#   include "setRootCase.H"
#   include "createTime.H"

    // Create beam model
    Foam::autoPtr<Foam::beamModel> beam =
        Foam::beamModel::New(runTime, Foam::dynamicFvMesh::defaultRegion);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Foam::Info<< "\nCalculating beam mean line displacement and"
              << " cross-section rotation\n" << Foam::endl;

    // while (runTime.run())
    while (runTime.loop())
    {
        // runTime.setDeltaT(beam().deltaT());

        // runTime++;

        Foam::Info<< "\n\nTime = " << runTime.timeName() << Foam::nl
                  << Foam::endl;

        beam().evolve();

        beam().updateTotalFields();

        beam().writeFields();
        // runTime.write();

        Foam::Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                  << Foam::nl << Foam::endl;
    }

    Foam::Info<< "End\n" << Foam::endl;

    return 0;
}


// ************************************************************************* //
