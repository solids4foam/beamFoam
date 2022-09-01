/*---------------------------------------------------------------------------*\
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

\*---------------------------------------------------------------------------*/

#include "conicalPulley.H"
#include "cylindricalCS.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// const dataType Foam::conicalPulley::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conicalPulley::conicalPulley()
:
    origin_(vector::zero),
    originShift_(vector::zero),
    temporalOriginDisplacement_(vector::zero),
    axis_(vector(0, 0, 1)),
    polarAxis_(vector(1, 0, 0)),
    rMin_(0),
    angle_(0),
    width_(0),
    rotDir_(1),
    meanBearingDiameter_(0),
    bearingFrictionCoeff_(0),
    frictionless_(false),
    frictionlessFullContact_(false),
    negConeApex_(vector::zero),
    posConeApex_(vector::zero),
    originDispSeries_(),
    stlFileName_("none"),
    stlModel_(),
    circumferentialSpeed_(0),
    angularVelocity_(0)
{}

Foam::conicalPulley::conicalPulley
(
    const vector& origin,
    const vector& originShift,
    const vector& axis,
    const vector& polarAxis,
    const scalar& rMin,
    const scalar& angle,
    const scalar& width,
    const label& rotDir
)
:
    origin_(origin),
    originShift_(originShift),
    temporalOriginDisplacement_(vector::zero),
    axis_(axis),
    polarAxis_(polarAxis),
    rMin_(rMin),
    angle_(angle),
    width_(width),
    rotDir_(rotDir),
    meanBearingDiameter_(0),
    bearingFrictionCoeff_(0),
    frictionless_(false),
    frictionlessFullContact_(false),
    negConeApex_(vector::zero),
    posConeApex_(vector::zero),
    originDispSeries_(),
    stlFileName_("none"),
    stlModel_(),
    circumferentialSpeed_(0),
    angularVelocity_(0)
{
    negConeApex_ = axis_*(rMin_*tan(angle_*M_PI/180));
    posConeApex_ = -axis_*(rMin_*tan(angle_*M_PI/180));
}

Foam::conicalPulley::conicalPulley(const dictionary& dict)
:
    origin_(dict.lookup("origin")),
    originShift_(dict.lookupOrDefault<vector>("originShift", vector::zero)),
    temporalOriginDisplacement_(vector::zero),
    axis_(dict.lookup("axis")),
    polarAxis_(dict.lookup("polarAxis")),
    rMin_(readScalar(dict.lookup("rMin"))),
    angle_(readScalar(dict.lookup("angle"))),
    width_(readScalar(dict.lookup("width"))),
    rotDir_(dict.lookupOrDefault<label>("rotDir", 1)),
    meanBearingDiameter_(dict.lookupOrDefault<scalar>("meanBearingDiameter", 0)),
    bearingFrictionCoeff_(dict.lookupOrDefault<scalar>("bearingFrictionCoeff", 0)),
    frictionless_(dict.lookupOrDefault<bool>("frictionless", false)),    
    frictionlessFullContact_(dict.lookupOrDefault<bool>("frictionlessFullContact", false)),    
    negConeApex_(vector::zero),
    posConeApex_(vector::zero),
    originDispSeries_(),
    stlFileName_(dict.lookupOrDefault<fileName>("stlFileName", "none")),
    stlModel_(),
    circumferentialSpeed_(dict.lookupOrDefault<scalar>("circumferentialSpeed", 0)),
    angularVelocity_(0)
{
    negConeApex_ = axis_*(rMin_*tan(angle_*M_PI/180));
    posConeApex_ = -axis_*(rMin_*tan(angle_*M_PI/180));

    // Check if origin displacement is time-varying
    if (dict.found("originDisplacementSeries"))
    {
        Info<< "Origin displacement is time-varying" << endl;
        originDispSeries_ =
            interpolationTable<vector>(dict.subDict("originDisplacementSeries"));
    }

    // Read stl model
    if (stlFileName_ != fileName("none"))
    {
        stlModel_ = triSurface(stlFileName_);
    }
}

// Foam::conicalPulley::conicalPulley(const conicalPulley&)
// :
//     data_()
// {}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Foam::autoPtr<Foam::conicalPulley> Foam::conicalPulley::New()
// {
//     return autoPtr<conicalPulley>(new conicalPulley);
// }

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conicalPulley::~conicalPulley()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::conicalPulley::position(const vector& cylCoord) const
{
    cylindricalCS ccs
    (
        "CCS",
        origin(),
        axis(),
        polarAxis(),
        false
    );

    vector pGlobalCartesian = ccs.globalPosition(cylCoord);

    // vector pLocalCartesian = vector::zero;

    // pLocalCartesian.z() = cylCoord.z();
    // pLocalCartesian.x() = cylCoord.x()*cos(cylCoord.y());
    // pLocalCartesian.y() = cylCoord.x()*sin(cylCoord.y());

    // // To do: transform into global cartesion coordinates
    // vector k(0, 0, 1);
    // vector i(1, 0, 0);
    // if
    // (
    //     (mag((k & axis_) - 1) > SMALL)
    //  || (mag((i & polarAxis_) - 1) > SMALL)
    // )
    // {
    //     FatalErrorIn("conicalPulley::position(...)")
    //         << "Locall frame of the pulley is not aligned with the global one."
    //         << "Transformation is required."
    //         << abort(FatalError);
    // }

    // // Global cartesian
    // pLocalCartesian += origin();

    return pGlobalCartesian;
}

Foam::vector Foam::conicalPulley::nearestPointNeg(const vector& p) const
{
    cylindricalCS ccs
    (
        "CCS",
        origin(),
        axis(),
        polarAxis(),
        false
    );
    
    vector pCyl = ccs.localPosition(p);

    // vector pCyl = vector::zero;

    // vector R = p - origin();
    // vector r = R - axis_*(axis_ & R);
    // vector rHat = r/mag(r);

    // scalar x = (polarAxis_ & rHat);
    // scalar y = ((axis_ ^ polarAxis_) & rHat);

    // // Radius
    // pCyl.x() = mag(r);
    // // Angle
    // pCyl.y() = atan2(y, x);
    // // Coordinate along axis
    // pCyl.z() = (R & axis_);

    vector nearest = vector::zero;

    if (pCyl.z() < -width()/2)
    {
        nearest.x() = maxRadius();
        nearest.y() = pCyl.y();
        nearest.z() = -width()/2;

        return nearest;
    }
    else if (pCyl.z() > width()/2)
    {
        nearest.x() = maxRadius();
        nearest.y() = pCyl.y();
        nearest.z() = width()/2;

        return nearest;
    }

    scalar rCone = rMin_ - pCyl.z()*tan(angle_*M_PI/180);
    // scala r d = (pCyl.x()-rCone)*cos(angle*M_PI/180);
    scalar c = (pCyl.x() - rCone)*sin(angle_*M_PI/180);
    scalar dz = c*cos(angle_*M_PI/180);
    scalar dr = c*sin(angle_*M_PI/180);

    nearest.x() = rCone + dr;
    if (nearest.x() < rMin_)
    {
        nearest.x() = rMin_;
    }
    
    nearest.y() = pCyl.y();

    nearest.z() = pCyl.z() - dz;
    if (nearest.z() > 0)
    {
        nearest.z() = 0;
    }
    // else if (nearest.z() < -width())
    // {
    //     nearest.z() = -width();
    // }
    
    return nearest;
}

Foam::vector Foam::conicalPulley::nearestPointPos(const vector& p) const
{
    cylindricalCS ccs
    (
        "CCS",
        origin(),
        axis(),
        polarAxis(),
        false
    );
    
    vector pCyl = ccs.localPosition(p);
    
    // vector pCyl = vector::zero;

    // vector R = p - origin();
    // vector r = R - axis_*(axis_ & R);
    // vector rHat = r/mag(r);

    // scalar x = (polarAxis_ & rHat);
    // scalar y = ((axis_ ^ polarAxis_) & rHat);

    // // Radius
    // pCyl.x() = mag(r);
    // // Angle
    // pCyl.y() = atan2(y, x);
    // // Coordinate along axis
    // pCyl.z() = (R & axis_);

    vector nearest = vector::zero;

    if (pCyl.z() > width()/2)
    {
        nearest.x() = maxRadius();
        nearest.y() = pCyl.y();
        nearest.z() = width()/2;

        return nearest;
    }
    else if (pCyl.z() < -width()/2)
    {
        nearest.x() = maxRadius();
        nearest.y() = pCyl.y();
        nearest.z() = -width()/2;

        return nearest;
    }
    
    scalar rCone = rMin_ + pCyl.z()*tan(angle_*M_PI/180); 
    // scalar d = (pCyl.x()-rCone)*cos(angle*M_PI/180);
    scalar c = (pCyl.x()-rCone)*sin(angle_*M_PI/180);
    scalar dz = c*cos(angle_*M_PI/180);
    scalar dr = c*sin(angle_*M_PI/180);

    nearest.x() = rCone + dr;
    if (nearest.x() < rMin_)
    {
        nearest.x() = rMin_;
    }
    
    nearest.y() = pCyl.y();

    nearest.z() = pCyl.z() + dz;
    if (nearest.z() < 0)
    {
        nearest.z() = 0;
    }

    return nearest;
}


void Foam::conicalPulley::move(const scalar time)
{
    if (originDispSeries_.size())
    {
        temporalOriginDisplacement_ = originDispSeries_(time);
    }
    else
    {
        temporalOriginDisplacement_ = vector::zero;
    }
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// void Foam::conicalPulley::operator=(const conicalPulley& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("Foam::conicalPulley::operator=(const Foam::conicalPulley&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// }

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
