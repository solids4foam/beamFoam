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

//- Monolithic coupling of a body to the attachment cell of the beam. Rows
//  and columns are Newton increments: body translation (force rows) and,
//  when rotation is active, body rotation (moment rows). Without rotation the
//  body's rotation rows are identity with zero right-hand side.
struct RigidBodyMonolithicCoupling
{
    bool active;

    bool rotationActive;

    //- Beam cell next to the attachment face
    label attachmentCell;

    //- Beam force rows of the attachment cell, body translation columns
    tensor beamWRowCoeff;

    //- Beam moment rows of the attachment cell, body translation columns
    tensor beamThetaRowCoeff;

    //- Body translation rows, body translation columns
    tensor bodyCoeff;

    //- Body translation rows, attachment cell displacement columns
    tensor bodyWCoeff;

    //- Body translation rows, attachment cell rotation columns
    tensor bodyThetaCoeff;

    //- Body translation right-hand side (minus the body residual)
    vector bodySource;

    //- Beam force and moment rows of the attachment cell, body rotation
    //  columns
    tensor beamWRowRotCoeff;
    tensor beamThetaRowRotCoeff;

    //- Body translation rows, body rotation columns
    tensor bodyTransRotCoeff;

    //- Body rotation rows: translation, rotation, attachment cell
    //  displacement and rotation columns
    tensor bodyRotTransCoeff;
    tensor bodyRotCoeff;
    tensor bodyRotWCoeff;
    tensor bodyRotThetaCoeff;

    //- Body rotation right-hand side (minus the rotational residual)
    vector bodyRotSource;

    RigidBodyMonolithicCoupling()
    :
        active(false),
        rotationActive(false),
        attachmentCell(-1),
        beamWRowCoeff(tensor::zero),
        beamThetaRowCoeff(tensor::zero),
        bodyCoeff(tensor::zero),
        bodyWCoeff(tensor::zero),
        bodyThetaCoeff(tensor::zero),
        bodySource(vector::zero),
        beamWRowRotCoeff(tensor::zero),
        beamThetaRowRotCoeff(tensor::zero),
        bodyTransRotCoeff(tensor::zero),
        bodyRotTransCoeff(tensor::zero),
        bodyRotCoeff(tensor::zero),
        bodyRotWCoeff(tensor::zero),
        bodyRotThetaCoeff(tensor::zero),
        bodyRotSource(vector::zero)
    {}
};

//- Accepted state of a rigidBodyEnd solved by beamFoam, for a caller such as
//  a sixDoF motion solver. Global axes; the body turns by orientation from
//  body to global axes
struct RigidBodyEndState
{
    vector displacement;
    tensor orientation;
    vector velocity;
    vector acceleration;
    vector angularVelocity;
    vector angularAcceleration;
    vector beamForce;
    vector beamTorque;
    label iterations;

    RigidBodyEndState()
    :
        displacement(vector::zero),
        orientation(tensor::I),
        velocity(vector::zero),
        acceleration(vector::zero),
        angularVelocity(vector::zero),
        angularAcceleration(vector::zero),
        beamForce(vector::zero),
        beamTorque(vector::zero),
        iterations(0)
    {}
};

}

#endif
