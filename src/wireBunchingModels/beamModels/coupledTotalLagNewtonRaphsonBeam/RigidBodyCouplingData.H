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

Description
    Data passed between coupledTotalLagNewtonRaphsonBeam and its rigidBodyEnd:
    RigidBodyMonolithicCoupling holds the body rows and columns of the
    monolithic beam + rigid body Newton system, and RigidBodyEndState the
    accepted body state for a caller such as a sixDoF motion solver.

Author
    Colm Mc Alister, UCD.

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
