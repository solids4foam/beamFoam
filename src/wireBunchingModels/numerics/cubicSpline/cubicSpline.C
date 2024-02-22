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

#include "error.H"
#include "cubicSpline.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cubicSpline::calcParam()
{
    param_.setSize(points_.size());

    if (param_.size())
    {
        param_[0] = 0.0;

        for (label i=1; i < param_.size(); i++)
        {
            param_[i] = param_[i-1] + mag(points_[i] - points_[i-1]);
        }

        lineLength_ = param_.last();

        // // normalize on the interval 0-1
        // for (label i=1; i < param_.size() - 1; i++)
        // {
        //     param_[i] /= lineLength_;
        // }
        // param_.last() = 1.0;
    }
    else
    {
        lineLength_ = 0.0;
    }
}

void Foam::cubicSpline::calcCubicSpline
(
    const bcType startBC,
    const vector startDerivative,
    const bcType endBC,
    const vector endDerivative
)
{
    scalarField a(points_.size(), 0);
    scalarField b(points_.size(), 0);
    scalarField c(points_.size(), 0);
    vectorField d(points_.size(), vector::zero);

    // vector D0(1, 0, 0);

    DD_[0] = vector::zero;
    DD_[nSegments()] = vector::zero;

    D_[0] = startDerivative;
    if (startBC == CLAMPED_CALC)
    {
        label i=1;
        scalar hp = mag(points_[i+1] - points_[i]);
        vector dRp = points_[i+1] - points_[i];
        scalar hm = mag(points_[i-1] - points_[i]);
        vector dRm = points_[i-1] - points_[i];
        vector c = (hm*dRp+hp*dRm)/(hp*hm*(hm+hp));
        vector b = dRp/hp - (hm*dRp+hp*dRm)/(hm*(hm+hp));

        D_[0] = b - 2*c*hm;
    }

    D_[nSegments()] = endDerivative;
    if (endBC == CLAMPED_CALC)
    {
        label i=nSegments()-1;
        scalar hp = mag(points_[i+1] - points_[i]);
        vector dRp = points_[i+1] - points_[i];
        scalar hm = mag(points_[i-1] - points_[i]);
        vector dRm = points_[i-1] - points_[i];
        vector c = (hm*dRp+hp*dRm)/(hp*hm*(hm+hp));
        vector b = dRp/hp - (hm*dRp+hp*dRm)/(hm*(hm+hp));

        D_[nSegments()] = b + 2*c*hp;
    }


    for (label i=1; i<(a.size()-1); i++)
    {
        a[i] = mag(points_[i] - points_[i-1]);
        b[i] = 2*mag(points_[i] - points_[i-1])
          + 2*mag(points_[i+1] - points_[i]);
        c[i] = mag(points_[i+1] - points_[i]);

        d[i] = 6*(points_[i+1]-points_[i])/mag(points_[i+1]-points_[i])
          - 6*(points_[i]-points_[i-1])/mag(points_[i]-points_[i-1]);

        if ( (i==1) && ( (startBC == CLAMPED) || (startBC == CLAMPED_CALC) ) )
        //ZT, Bug fix - CLAMPED -> CLAMPED_CALC
        {
            b[i] -= 0.5*mag(points_[i] - points_[i-1]);
            d[i] += 3*D_[0]
              - 3*(points_[i]-points_[i-1])/mag(points_[i]-points_[i-1]);
        }

        if
        (
            (i == (a.size()-2))
         && ( (endBC == CLAMPED) || (endBC == CLAMPED_CALC) )
        //ZT, Bug fix - CLAMPED -> CLAMPED_CALC
        )
        {
            b[i] -= 0.5*mag(points_[i+1] - points_[i]);

            d[i] += 3*(points_[i+1]-points_[i])/mag(points_[i+1]-points_[i])
              - 3*D_[nSegments()];
        }
    }

    // Solve three-diagonal system for second derivative
    for (label i=1; i<(a.size()-2); i++)
    {
        if (i==1)
        {
            c[i] = c[i]/b[i];
        }
        else
        {
            c[i] = c[i]/(b[i]-a[i]*c[i-1]);
        }
    }

    for (label i=1; i<(a.size()-1); i++)
    {
        if (i==1)
        {
            d[i] = d[i]/b[i];
        }
        else
        {
            d[i] = (d[i]-a[i]*d[i-1])/(b[i]-a[i]*c[i-1]);
        }
    }

    // Back-substitution
    for (label i=(a.size()-2); i>0; i--)
    {
        if (i==(a.size()-2))
        {
            DD_[i] = d[i];
        }
        else
        {
            DD_[i] = d[i] - c[i]*DD_[i+1];
        }
    }

    if ( (startBC == CLAMPED) || (startBC == CLAMPED_CALC) )
    {
        DD_[0] = 3*(points_[1]-points_[0])/magSqr(points_[1]-points_[0])
          - 0.5*DD_[1] - 3*D_[0]/mag(points_[1]-points_[0]);
    }

    if ( (endBC == CLAMPED) || (endBC == CLAMPED_CALC) )
    {
        label i = DD_.size()-1;
        DD_[i] = -3*(points_[i]-points_[i-1])/magSqr(points_[i]-points_[i-1])
          - 0.5*DD_[i-1] + 3*D_[nSegments()]/mag(points_[i]-points_[i-1]);
    }

    // First derivative dRdt (!= dRdS)
    for (label i=0; i<nSegments(); i++)
    {
        scalar hi = mag(points_[i+1] - points_[i]);
        D_[i] = -DD_[i]*hi/2 + (points_[i+1] - points_[i])/hi
          - (DD_[i+1] - DD_[i])*hi/6;
    }
    label i = nSegments();
    scalar hi = mag(points_[i] - points_[i-1]);
    D_[i] = DD_[i]*hi/2 + (points_[i] - points_[i-1])/hi
      - (DD_[i] - DD_[i-1])*hi/6;
}

void Foam::cubicSpline::calcMidPoints() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (midPointsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::cubicSpline::calcMidPoints() const"
        )
            << "Segment mid points already exist"
            << abort(FatalError);
    }

    midPointsPtr_ = new vectorField(nSegments(), vector::zero);
    vectorField& midPoints = *midPointsPtr_;

    forAll(midPoints, sI)
    {
        midPoints[sI] = position(sI, 0.5);
    }
}

void Foam::cubicSpline::calcMidPointDerivatives() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (midPointDerivativesPtr_)
    {
        FatalErrorIn
        (
            "void Foam::cubicSpline::calcMidPointDerivatives() const"
        )
            << "Segment mid point derivatives already exist"
            << abort(FatalError);
    }

    midPointDerivativesPtr_ = new vectorField(nSegments(), vector::zero);
    vectorField& midPointDerivatives = *midPointDerivativesPtr_;

    forAll(midPointDerivatives, sI)
    {
        midPointDerivatives[sI] = firstDerivative(sI, 0.5);
    }

    midPointDerivatives /= mag(midPointDerivatives);
}

void Foam::cubicSpline::calcDRdS() const
{
    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (dRdSPtr_)
    {
        FatalErrorIn
        (
            "void Foam::cubicSpline::calcDRdS() const"
        )
            << "Segment end point derivatives already exist"
            << abort(FatalError);
    }

    dRdSPtr_ = new vectorField(D_);
    vectorField& dRdS = *dRdSPtr_;

    dRdS /= mag(dRdS);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cubicSpline::cubicSpline
(
    const pointField& ps,
    const bcType startBC,
    const vector startDerivative,
    const bcType endBC,
    const vector endDerivative
)
:
    points_(ps),
    midPointsPtr_(NULL),
    midPointDerivativesPtr_(NULL),
    dRdSPtr_(NULL),
    lineLength_(0.0),
    param_(0),
    D_(ps.size(), vector::zero),
    DD_(ps.size(), vector::zero)
{
    calcParam();
    calcCubicSpline(startBC, startDerivative, endBC, endDerivative);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::cubicSpline::points() const
{
    return points_;
}


Foam::label Foam::cubicSpline::nSegments() const
{
    return points_.size()-1;
}


Foam::label Foam::cubicSpline::localParameter(scalar& lambda) const
{
    // check endpoints
    if (lambda < SMALL)
    {
        lambda = 0;
        return 0;
    }
    else if (lambda > lineLength_ - SMALL)
    {
        lambda = 1;
        return nSegments();
    }

    // search table of cumulative distances to find which line-segment
    // we are on. Check the upper bound.

    label segmentI = 1;
    while (param_[segmentI] < lambda)
    {
        segmentI++;
    }
    segmentI--;   // we want the corresponding lower bound

    // the local parameter [0-1] on this line segment
    lambda =
    (
        ( lambda - param_[segmentI] )
      / ( param_[segmentI+1] - param_[segmentI] )
    );

    return segmentI;
}


Foam::point Foam::cubicSpline::position(const scalar mu) const
{
    // check endpoints
    if (mu < SMALL)
    {
        return points_.first();
    }
    else if (mu > lineLength_ - SMALL)
    {
        return points_.last();
    }

    scalar lambda = mu;
    label segment = localParameter(lambda);
    return position(segment, lambda);
}


Foam::point Foam::cubicSpline::position
(
    const label segment,
    const scalar mu
) const
{
    // out-of-bounds
    if (segment < 0)
    {
        return points_.first();
    }
    else if (segment > nSegments())
    {
        return points_.last();
    }

    const point& p0 = points()[segment];
    const point& p1 = points()[segment+1];

    // special cases - no calculation needed
    if (mu <= 0.0)
    {
        return p0;
    }
    else if (mu >= 1.0)
    {
        return p1;
    }
    else
    {
        // cubic spline interpolation

        scalar h = param_[segment+1] - param_[segment];
        scalar t = param_[segment] + mu*h;
        scalar h0 = t - param_[segment];
        scalar h1 = param_[segment+1] - t;

        vector curR =
            DD_[segment+1]*Foam::pow(h0, 3)/(6*h)
          + DD_[segment]*Foam::pow(h1, 3)/(6*h)
          + (points_[segment+1]/h - DD_[segment+1]*h/6)*h0
          + (points_[segment]/h - DD_[segment]*h/6)*h1;

        // // Linear interpolation
        // curR = points_[segment]
        //   + (
        //         (points_[segment+1]-points_[segment])
        //        /(param_[segment+1]-param_[segment])
        //     )
        //    *(t-param_[segment]);

        return curR;
    }
}


Foam::point Foam::cubicSpline::firstDerivative(const scalar mu) const
{
    // check endpoints
    if (mu < SMALL)
    {
        return D_.first();
    }
    else if (mu > lineLength_ - SMALL)
    {
        return D_.last();
    }

    scalar lambda = mu;
    label segment = localParameter(lambda);
    return firstDerivative(segment, lambda);
}

Foam::point Foam::cubicSpline::firstDerivative
(
    const label segment,
    const scalar mu
) const
{
    // out-of-bounds
    if (segment < 0)
    {
        return D_.first();
    }
    else if (segment > nSegments())
    {
        return D_.last();
    }

    const point& D0 = D_[segment];
    const point& D1 = D_[segment+1];

    // special cases - no calculation needed
    if (mu <= 0.0)
    {
        return D0;
    }
    else if (mu >= 1.0)
    {
        return D1;
    }
    else
    {
        // cubic spline interpolation

        scalar h = param_[segment+1] - param_[segment];
        scalar t = param_[segment] + mu*h;
        scalar h0 = t - param_[segment];
        scalar h1 = param_[segment+1] - t;

        vector curD =
            DD_[segment+1]*Foam::pow(h0, 2)/(2*h)
          - DD_[segment]*Foam::pow(h1, 2)/(2*h)
          + (points_[segment+1] - points_[segment])/h
          + (DD_[segment+1] - DD_[segment]);

        return curD;
    }
}

Foam::scalar Foam::cubicSpline::length() const
{
    return lineLength_;
}

const Foam::vectorField& Foam::cubicSpline::midPoints() const
{
    if (!midPointsPtr_)
    {
        calcMidPoints();
    }

    return *midPointsPtr_;
}

const Foam::vectorField& Foam::cubicSpline::dRdS() const
{
    if (!dRdSPtr_)
    {
        calcDRdS();
    }

    return *dRdSPtr_;
}

const Foam::vectorField& Foam::cubicSpline::midPointDerivatives() const
{
    if (!midPointDerivativesPtr_)
    {
        calcMidPointDerivatives();
    }

    return *midPointDerivativesPtr_;
}

Foam::labelScalar Foam::cubicSpline::nearestPoint
// Foam::vector Foam::cubicSpline::nearestPoint
(
    const label segI,
    const vector& p
) const
{
    vector np(0, 0, 0);

    scalar mu0 = 0;
    scalar f0 = (D_[segI] & (p - points_[segI]));

    scalar mu1 = 1;
    scalar f1 = (D_[segI+1] & (p - points_[segI+1]));

    scalar mu = 0;

    if ((f0*f1) > SMALL)
    {
        FatalErrorIn
        (
            "Foam::cubicSpline::nearestPoint(...) const"
        )
            << "Nearest point is out of segment"
            << abort(FatalError);
    }
    else if (mag(f0) < SMALL)
    {
        np = points_[segI];
    }
    else if (mag(f1) < SMALL)
    {
        np = points_[segI+1];
    }
    else // find root
    {
        scalar err = GREAT;
        do
        {
            scalar df1 = (f1-f0)/(mu1-mu0);
            mu = mu1 - (f1/df1);
            scalar f =
            (
                firstDerivative(segI, mu)
              & (p - position(segI, mu))
            );

            if (f0*f>0)
            {
                err = mag(mu-mu0);
                mu0 = mu;
                f0 = f;
            }
            else
            {
                err = mag(mu-mu1);
                mu1 = mu;
                f1 = f;
            }
        }
        while(err<1e-3);

        np = position(segI, mu);
    }

    return labelScalar(segI, mu);
}

Foam::tmp<Foam::scalarField> Foam::cubicSpline::segLengths() const
{
    tmp<scalarField> tSegLengths
    (
        new scalarField(nSegments(), 0)
    );

    forAll(tSegLengths(), segI)
    {
        tSegLengths.ref()[segI] = mag(points_[segI+1] - points_[segI]);
    }

    return tSegLengths;
}

Foam::tmp<Foam::scalarField> Foam::cubicSpline::midPointParameters() const
{
    tmp<scalarField> tMidParameters
    (
        new scalarField(nSegments(), 0)
    );

    const scalarField segLen (segLengths());

    forAll(tMidParameters(), segI)
    {
        tMidParameters.ref()[segI] = param_[segI] + segLen[segI]/2;
    }

    return tMidParameters;
}


Foam::tmp<Foam::vectorField> Foam::cubicSpline::segmentToPointInterpolate
(
    const vectorField& midPhi
) const
{
    tmp<vectorField> tResult
    (
        new vectorField(points_.size(), vector::zero)
    );
    vectorField& result = tResult.ref();

    label n = nSegments() + 2;

    scalarField a(n, 0);
    scalarField b(n, 0);
    scalarField c(n, 0);
    vectorField d(n, vector::zero);

    vectorField D(n, vector::zero);
    vectorField DD(n, vector::zero);

    vectorField points = midPhi;

    scalarField midPointParam (midPointParameters());

    const bcType startBC = CLAMPED;
    // const bcType startBC = CLAMPED_CALC;
    const vector startDerivative = vector::zero;
    const bcType endBC = CLAMPED;
    // const bcType endBC = CLAMPED_CALC;
    const vector endDerivative = vector::zero;

    DD[0] = vector::zero;
    DD[n-1] = vector::zero;

    D[0] = startDerivative;
    if (startBC == CLAMPED_CALC)
    {
        label i=1;
        scalar hp = midPointParam[i] - midPointParam[i-1];
        vector dRp = points[i+1] - points[i];
        scalar hm = midPointParam[i-1];
        vector dRm = points[i-1] - points[i];
        vector c = (hm*dRp+hp*dRm)/(hp*hm*(hm+hp));
        vector b = dRp/hp - (hm*dRp+hp*dRm)/(hm*(hm+hp));

        D[0] = b - 2*c*hm;
    }

    D[n-1] = endDerivative;
    if (endBC == CLAMPED_CALC)
    {
        label i=n-2;
        scalar hp = param_.last() - midPointParam[i-1];
        vector dRp = points[i+1] - points[i];
        scalar hm = midPointParam[i-1] - midPointParam[i-2];
        vector dRm = points[i-1] - points[i];
        vector c = (hm*dRp+hp*dRm)/(hp*hm*(hm+hp));
        vector b = dRp/hp - (hm*dRp+hp*dRm)/(hm*(hm+hp));

        D[n-1] = b + 2*c*hp;
    }

    for (label i=1; i<(a.size()-1); i++)
    {
        if (i==1)
        {
            a[i] = midPointParam[i-1] - 0;
            c[i] = midPointParam[i] - midPointParam[i-1];
        }
        else if (i==(a.size()-2))
        {
            a[i] = midPointParam[i-1] - midPointParam[i-2];
            c[i] = param_.last() - midPointParam[i-1];
        }
        else
        {
            a[i] = midPointParam[i-1] - midPointParam[i-2];
            c[i] = midPointParam[i] - midPointParam[i-1];
        }

        b[i] = 2*a[i] + 2*c[i];
        d[i] = 6*(points[i+1]-points[i])/c[i]
          - 6*(points[i]-points[i-1])/a[i];

        if
        (
            (i==1)
         && (
               (startBC == CLAMPED)
            || (startBC == CLAMPED_CALC)
            )
        )
        {
            b[i] -= 0.5*midPointParam[i-1];
            d[i] += 3*D[0] - 3*(points[i]-points[i-1])/midPointParam[i-1];
        }

        if
        (
            (i == (a.size()-2))
         && (
                (endBC == CLAMPED)
             || (endBC == CLAMPED_CALC)  // Bug fix?
            )
        )
        {
            b[i] -= 0.5*(param_.last() - midPointParam[i-1]);
            d[i] += 3*(points[i+1]-points[i])
               /(param_.last() - midPointParam[i-1]) - 3*D[n-1];
        }
    }

    // Solve three-diagonal system for second derivative
    for (label i=1; i<(a.size()-2); i++)
    {
        if (i==1)
        {
            c[i] = c[i]/b[i];
        }
        else
        {
            c[i] = c[i]/(b[i]-a[i]*c[i-1]);
        }
    }

    for (label i=1; i<(a.size()-1); i++)
    {
        if (i==1)
        {
            d[i] = d[i]/b[i];
        }
        else
        {
            d[i] = (d[i]-a[i]*d[i-1])/(b[i]-a[i]*c[i-1]);
        }
    }

    // Back-substitution
    for (label i=(a.size()-2); i>0; i--)
    {
        if (i==(a.size()-2))
        {
            DD[i] = d[i];
        }
        else
        {
            DD[i] = d[i] - c[i]*DD[i+1];
        }
    }

    if ( (startBC == CLAMPED) || (startBC == CLAMPED_CALC) )
    {
        DD[0] = 3*(points[1]-points[0])/sqr(midPointParam[0])
          - 0.5*DD[1] - 3*D[0]/midPointParam[0];
    }

    if ( (endBC == CLAMPED) || (endBC == CLAMPED_CALC) )
    {
        label i = DD.size()-1;
        DD[i] = -3*(points[i]-points[i-1])
           /sqr(param_.last() - midPointParam[i-1])
          - 0.5*DD[i-1] + 3*D[i]/sqr(param_.last() - midPointParam[i-1]);
    }

    // Interpolate
    result[0] = points[0];
    result[result.size()-1] = points[points.size()-1];
    for (label i=1; i<(points.size()-1); i++)
    {
        scalar h = midPointParam[i] - midPointParam[i-1];
        scalar t = midPointParam[i-1] + 0.5*h;
        scalar h0 = t - midPointParam[i-1];
        scalar h1 = midPointParam[i] - t;

        result[i] = DD_[i+1]*Foam::pow(h0, 3)/(6*h)
          + DD_[i]*Foam::pow(h1, 3)/(6*h)
          + (points[i+1]/h - DD_[i+1]*h/6)*h0
          + (points[i]/h - DD_[i]*h/6)*h1;
    }

    return tResult;
}


// ************************************************************************* //
