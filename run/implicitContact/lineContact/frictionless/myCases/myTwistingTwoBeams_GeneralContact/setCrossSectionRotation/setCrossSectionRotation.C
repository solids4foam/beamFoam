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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "setCrossSectionRotation.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"
#include "axialForceTransverseDisplacementFvPatchVectorField.H"
#include "pseudoVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setCrossSectionRotation, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setCrossSectionRotation,
        dictionary
    );
}


Foam::tensor R(const Foam::tensor& T0, const Foam::tensor& T1);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setCrossSectionRotation::setBC()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    scalar tEnd = time_.endTime().value();

    // vectorField curPosition =
    //     mesh.C().boundaryField()[startPatchIndex_];

    // vector i(1, 0, 0);
    // scalarField curRadius = mag((I-i*i) & curPosition);
    // Info << regionName_ << ", curRadius: " << curRadius << endl;
        
    // if (mesh.foundObject<volVectorField>("DTheta"))
    // {
    // }
    // else
    {
        volVectorField& Theta =
            const_cast<volVectorField&>
            (
                mesh.lookupObject<volVectorField>("Theta")
            );

        // scalar oldTheta =
        //     thetaEnd_*(time_.value()/tEnd);
        
        scalar theta =
            thetaEnd_*((time_.value() + time_.deltaT().value())/tEnd);

        scalar r = 0.095;
        scalar h = 5.049647;
        scalar hPrime = 2*M_PI*h/theta;

        // tensor T0
        // (
        //     1,  0,  0,
        //     0, -1,  0,
        //     0,  0, -1
        // );
        
        vector t0(1, 0, 0);

        // Start points
        {
            scalar theta = 0;

            // First beam
            {
                if (mesh.moving()) // Updated Lagrangian formulation
                {
                    scalar oldTheta = thetaEnd_*time_.value()/tEnd;
                    scalar hPrime = 2*M_PI*h/oldTheta;

                    t0 =
                       vector
                       (
                           hPrime,
                          -r*::sin(theta),
                           r*::cos(theta)
                       );
                    t0 /= mag(t0);
                }

                vector t
                (
                    hPrime,
                   -r*::sin(theta),
                    r*::cos(theta)
                );
                t /= mag(t);

                scalar angleScalar = ::acos(t0 & t);
                vector axis = (t0 ^ t);
                axis /= mag(axis);
                vector angle = angleScalar*axis;
                
                Theta.boundaryField()[startPatchIndex_][0] == angle;
            }

            // Second beam
            if (true)
            {
                if (mesh.moving()) // Updated Lagrangian formulation
                {
                    scalar oldTheta = thetaEnd_*time_.value()/tEnd;
                    scalar hPrime = 2*M_PI*h/oldTheta;

                    t0 =
                       vector
                       (
                           hPrime,
                          -r*::sin(theta + M_PI),
                           r*::cos(theta + M_PI)
                       );
                    t0 /= mag(t0);
                }
                
                vector t
                (
                    hPrime,
                   -r*::sin(theta + M_PI),
                    r*::cos(theta + M_PI)
                );
                t /= mag(t);

                scalar angleScalar = ::acos(t0 & t);
                vector axis = (t0 ^ t);
                axis /= mag(axis);
                vector angle = angleScalar*axis;

                Theta.boundaryField()[startPatchIndex_][1] == angle;
            }
        }
        
        // End points
        {
            // First beam
            {
                if (mesh.moving()) // Updated Lagrangian formulation
                {
                    scalar oldTheta = thetaEnd_*time_.value()/tEnd;
                    scalar hPrime = 2*M_PI*h/oldTheta;

                    t0 =
                       vector
                       (
                           hPrime,
                          -r*::sin(oldTheta),
                           r*::cos(oldTheta)
                       );
                    t0 /= mag(t0);
                }
                
                vector t
                (
                    hPrime,
                   -r*::sin(theta),
                    r*::cos(theta)
                );
                t /= mag(t);

                scalar angleScalar = ::acos(t0 & t);
                vector axis = (t0 ^ t);
                axis /= mag(axis);
                vector angle = angleScalar*axis;

                Theta.boundaryField()[endPatchIndex_][0] == angle;  
            }
            
            // Second beam
            if (true)
            {
                if (mesh.moving()) // Updated Lagrangian formulation
                {
                    scalar oldTheta = thetaEnd_*time_.value()/tEnd;
                    scalar hPrime = 2*M_PI*h/oldTheta;

                    t0 =
                       vector
                       (
                           hPrime,
                          -r*::sin(oldTheta + M_PI),
                           r*::cos(oldTheta + M_PI)
                       );
                    t0 /= mag(t0);
                }
                
                vector t
                (
                    hPrime,
                   -r*::sin(theta + M_PI),
                    r*::cos(theta + M_PI)
                );
                t /= mag(t);

                scalar angleScalar = ::acos(t0 & t);
                vector axis = (t0 ^ t);
                axis /= mag(axis);
                vector angle = angleScalar*axis;
                
                Theta.boundaryField()[endPatchIndex_][1] == angle;
            }
        }
    }
    
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setCrossSectionRotation::setCrossSectionRotation
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    startPatchIndex_(-1),
    endPatchIndex_(-1),
    radius_(-1),
    // radius_(readScalar(dict.lookup("radius"))),
    thetaEnd_(readScalar(dict.lookup("angle"))*M_PI/180)
{
    Info << "Creating setCrossSectionRotation function object" << endl;
    
    if (Pstream::parRun())
    {
        FatalErrorIn("setCrossSectionRotation::setCrossSectionRotation(...)")
            << "setCrossSectionRotation objec function "
                << "is not implemented for parallel run"
                << abort(FatalError);
    }

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    word startPatchName(dict.lookup("startPatchName"));

    polyPatchID startPatch(startPatchName, mesh.boundaryMesh());
    
    if (!startPatch.active())
    {
        FatalErrorIn("setCrossSectionRotation::setCrossSectionRotation(...)")
          << "Patch name " << startPatchName << " not found."
          << abort(FatalError);
    }

    startPatchIndex_ = startPatch.index();

    word endPatchName(dict.lookup("endPatchName"));

    polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

    if (!endPatch.active())
    {
        FatalErrorIn("setCrossSectionRotation::setCrossSectionRotation(...)")
            << "Patch name " << endPatchName << " not found."
            << abort(FatalError);
    }

    endPatchIndex_ = endPatch.index();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setCrossSectionRotation::start()
{
    return setBC();
}

#if FOAMEXTEND > 40
bool Foam::setCrossSectionRotation::execute(const bool forceWrite)
#else
bool Foam::setCrossSectionRotation::execute()
#endif
{
    return setBC();
}

bool Foam::setCrossSectionRotation::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


inline Foam::tensor R
(
    const Foam::tensor& T0,
    const Foam::tensor& T1
)
{
    Foam::vector X0(T0.xx(), T0.yx(), T0.zx());
    Foam::vector Y0(T0.xy(), T0.yy(), T0.zy());
    Foam::vector Z0(T0.xz(), T0.yz(), T0.zz());

    Foam::vector X1(T1.xx(), T1.yx(), T1.zx());
    Foam::vector Y1(T1.xy(), T1.yy(), T1.zy());
    Foam::vector Z1(T1.xz(), T1.yz(), T1.zz());

    // tensor result
    // (
    //    (X0&X1), (X0&Y1), (X0&Z0), 
    //    (Y0&X1), (Y0&Y1), (Y0&Z0), 
    //    (Z0&X1), (Z0&Y1), (Z0&Z0)      
    // );
    Foam::tensor result
    (
       (X1&X0), (X1&Y0), (X1&Z0), 
       (Y1&X0), (Y1&Y0), (Y1&Z0), 
       (Z1&X0), (Z1&Y0), (Z1&Z0)      
    );

    result = (T1 & (T1 & result));
    
    return result;
}


// ************************************************************************* //
