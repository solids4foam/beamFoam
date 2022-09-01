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

#include "toroidalPulley.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// const dataType Foam::toroidalPulley::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::toroidalPulley::toroidalPulley()
:
    origin_(vector::zero),
    originShift_(vector::zero),
    temporalOriginDisplacement_(vector::zero),
    axis_(vector(0, 0, 1)),
    polarAxis_(vector(1, 0, 0)),
    rMin_(0),
    rToroidalSurf_(0),
    width_(0),
    rotDir_(1),
    originDispSeries_(),
    stlFileName_("none"),
    stlModel_()
{}

Foam::toroidalPulley::toroidalPulley
(
    const vector& origin,
    const vector& originShift,
    const vector& axis,
    const vector& polarAxis,
    const scalar& rMin,
    const scalar& rToroidalSurf,
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
    rToroidalSurf_(rToroidalSurf),
    width_(width),
    rotDir_(rotDir),
    meanBearingDiameter_(0),
    bearingFrictionCoeff_(0),
    originDispSeries_(),
    stlFileName_("none"),
    stlModel_()
{}

Foam::toroidalPulley::toroidalPulley(const dictionary& dict)
:
    origin_(dict.lookup("origin")),
    originShift_(dict.lookupOrDefault<vector>("originShift", vector::zero)),
    temporalOriginDisplacement_(vector::zero),
    axis_(dict.lookup("axis")),
    polarAxis_(dict.lookup("polarAxis")),
    rMin_(readScalar(dict.lookup("rMin"))),
    rToroidalSurf_(readScalar(dict.lookup("rToroidalSurf"))),
    width_(readScalar(dict.lookup("width"))),
    rotDir_(dict.lookupOrDefault<label>("rotDir", 1)),
    meanBearingDiameter_(dict.lookupOrDefault<scalar>("meanBearingDiameter", 0)),
    bearingFrictionCoeff_(dict.lookupOrDefault<scalar>("bearingFrictionCoeff", 0)),
    originDispSeries_(),
    stlFileName_(dict.lookupOrDefault<fileName>("stlFileName", "none")),
    stlModel_()
{
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

// Foam::toroidalPulley::toroidalPulley(const toroidalPulley&)
// :
//     data_()
// {}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Foam::autoPtr<Foam::toroidalPulley> Foam::toroidalPulley::New()
// {
//     return autoPtr<toroidalPulley>(new toroidalPulley);
// }

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::toroidalPulley::~toroidalPulley()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::vector Foam::toroidalPulley::position(const vector& cylCoord) const
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

    return pGlobalCartesian;
}

Foam::vector Foam::toroidalPulley::nearestPoint(const vector& p) const
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

    vector nearest = vector::zero;

    if (pCyl.z() < -width()/2)
    {
        Info << "xxxxxxxxx" << endl;
        nearest.x() = maxRadius();
        nearest.y() = pCyl.y();
        nearest.z() = -width()/2;

    }
    else if (pCyl.z() > width()/2)
    {
        Info << "yyyyyyyyy" << endl;
        nearest.x() = maxRadius();
        nearest.y() = pCyl.y();
        nearest.z() = width()/2;
        
        return nearest;
    }

    if (rToroidalSurf_/rMin_ > 1000)
    {
        // Cylindrical pulley
        nearest.x() = minRadius();
        nearest.y() = pCyl.y();
        nearest.z() = pCyl.z();
        
        return nearest;
    }

    Info << "Not cylindrical pulley" << endl;
    
    vector torSurfCentre
    (
        rMin_ + rToroidalSurf_,
        pCyl.y(),
        0 // z
    );

    vector dr =
        ccs.globalPosition(pCyl)
      - ccs.globalPosition(torSurfCentre);
        
    if (mag(dr) < SMALL)
    {
        nearest = vector(rMin_, pCyl.y(), 0);

        return nearest;
    }
    else
    {
        dr /= mag(dr);
        dr *= rToroidalSurf_;
    }
            
    if (pCyl.z() < 0)
    {
        scalar rTop =
            maxRadius()
          + (torSurfCentre.x() - maxRadius())
           *(width()/2 + pCyl.z())/(width()/2);

        if (pCyl.x() <= rTop)
        {
            nearest = 
                ccs.localPosition
                (
                    ccs.globalPosition(torSurfCentre) + dr
                );
        }
        else
        {
            nearest.x() = maxRadius();
            nearest.y() = pCyl.y();
            nearest.z() = -width()/2;
        }
    }
    else
    {
        scalar rTop =
            torSurfCentre.x()
          + (maxRadius() - torSurfCentre.x())*pCyl.z()
           /(width()/2);

        if (pCyl.x() <= rTop)
        {
            nearest =
                ccs.localPosition
                (
                    ccs.globalPosition(torSurfCentre) + dr
                );
        }
        else
        {
            nearest.x() = maxRadius();
            nearest.y() = pCyl.y();
            nearest.z() = width()/2;
        }
    }
        
    // if (nearest.x() < rMin_)
    // {
    //     nearest.x() = rMin_;
    // }
    
    // nearest.y() = pCyl.y();
    // nearest.z() = pCyl.z() - dz;

    
    return nearest;
}


Foam::vector Foam::toroidalPulley::axialTangent(const vector& cylCoord) const
{
    cylindricalCS ccs
    (
        "CCS",
        origin(),
        axis(),
        polarAxis(),
        false
    );

    vector T = vector::zero;

    vector torSurfCentre
    (
        rMin_ + rToroidalSurf_,
        cylCoord.y(),
        0 // z
    );

    vector dr =
        ccs.globalPosition(cylCoord)
      - ccs.globalPosition(torSurfCentre);

    dr /= mag(dr) + SMALL;
    dr = ccs.localPosition(dr);

    T = vector(dr.z(), cylCoord.y(), -dr.x());

    // Cylindrical pulley
    T = vector(0, 0, -1);
    
    T = ccs.globalPosition(T);
    T /= mag(T) + SMALL;

    return T;
}


Foam::scalar Foam::toroidalPulley::DrDz(const vector& cylCoord) const
{
    scalar drdz = 0;

    if (rToroidalSurf_/rMin_ <= 1000)
    {
        drdz =
            cylCoord.z()
           /::sqrt(sqr(rToroidalSurf_) - sqr(cylCoord.z()));
    }

    return drdz;
}


void Foam::toroidalPulley::move(const scalar time)
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

// void Foam::toroidalPulley::operator=(const toroidalPulley& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("Foam::toroidalPulley::operator=(const Foam::toroidalPulley&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// }

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
