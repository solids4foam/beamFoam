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

#ifndef RigidBodyCouplingData_H
#define RigidBodyCouplingData_H

#include "vector.H"
#include "tensor.H"

namespace Foam
{

struct RigidBodyState
{
    vector displacement;
    tensor orientation;
    vector velocity;
    vector angularMomentum;
    vector acceleration;
    vector torque;

    RigidBodyState()
    :
        displacement(0, 0, 0),
        orientation(tensor::I),
        velocity(0, 0, 0),
        angularMomentum(0, 0, 0),
        acceleration(0, 0, 0),
        torque(0, 0, 0)
    {}
};

struct RigidBodyStepData
{
    RigidBodyState previous;
    RigidBodyState current;
};

struct RigidBodySolution
{
    vector displacement;
    vector rotationCorrection;

    RigidBodySolution()
    :
        displacement(0, 0, 0),
        rotationCorrection(0, 0, 0)
    {}
};

}

#endif
