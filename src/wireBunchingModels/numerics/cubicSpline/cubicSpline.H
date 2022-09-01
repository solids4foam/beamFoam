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

Class
    Foam::cubicSpline

Description
    A series of straight line segments, which can also be interpreted as
    a series of control points for splines, etc.

SourceFiles
    cubicSpline.C

\*---------------------------------------------------------------------------*/

#ifndef cubicSpline_H
#define cubicSpline_H

#include "pointField.H"
#include "scalarList.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class cubicSpline Declaration
\*---------------------------------------------------------------------------*/

typedef Tuple2<label, scalar> labelScalar;
  
class cubicSpline
{
public:

    // Public enumerations
  
        //- Types of boundary conditions
        enum bcType
        {
            NATURAL,
            CLAMPED,
            CLAMPED_CALC
        };

private:
  
    // Private Member Functions

        //- Disallow default bitwise copy construct
        cubicSpline(const cubicSpline&);

        //- Disallow default bitwise assignment
        void operator=(const cubicSpline&);

protected:

    // Protected data

        //- The end points of each segments
        pointField points_;

        //- The mid points of segments
        mutable vectorField* midPointsPtr_;

        //- The segement mid point derivatives
        mutable vectorField* midPointDerivativesPtr_;

        //- The derivative at end points
        mutable vectorField* dRdSPtr_;
  
        //- The real line length
        scalar lineLength_;

        //- The rational (0-1) cumulative parameter value for each point
        scalarList param_;

        //- First derivatives
        pointField D_;
  
        //- Second derivatives
        pointField DD_;
  
    // Protected Member Functions

        //- Precalculate the rational cumulative parameter value
        //  and the line-length
        void calcParam();

        //- Calculate second derivative by solvine three-diagonal
        //  system of linear equations
        void calcCubicSpline
        (
            const bcType startBC = NATURAL,
            const vector startDerivative = vector(0,0,0),
            const bcType endBC = NATURAL,
            const vector endDerivative = vector(0,0,0)
        );
  
        //- Return the line segment and the local parameter [0..1]
        //  corresponding to the global lambda [0..1]
        label localParameter(scalar& lambda) const;

        //- Calculate segment mid points
        void calcMidPoints() const;
  
        //- Calculate segment mid points
        void calcMidPointDerivatives() const;
  
        //- Calculate segment end point derivatives
        void calcDRdS() const;

public:

    // Constructors

        //- Construct from components
        cubicSpline
        (
            const pointField& points,
            const bcType startBC = NATURAL,
            const vector startDerivative = vector(0, 0, 0),
            const bcType endBC = NATURAL,
            const vector endDerivative = vector(0, 0, 0)
        );

    // Member Functions

        //- Return const-access to the control-points
        const pointField& points() const;

        //- Return the number of line segments
        label nSegments() const;

        //- Return the point position corresponding to the curve parameter
        //  (cord length)
        point position(const scalar) const;

        //- Return the point position corresponding to the local parameter
        point position(const label segment, const scalar) const;

        //- Return the first derivative corresponding to the curve parameter
        //  (cord length)
        point firstDerivative(const scalar) const;

        //- Return the first derivative corresponding to the local parameter
        point firstDerivative(const label segment, const scalar) const;

        //- Return the length of the curve
        scalar length() const;
  
        //- Return midpoints of spline segments
        const vectorField& midPoints() const;
  
        //- Return segment midpoint derivatives
        const vectorField& midPointDerivatives() const;
  
        //- Return tangents of spline segments
        const vectorField& dRdS() const;
  
        //- Nearest segment point for given point
        labelScalar nearestPoint(const label segI, const vector& point) const;
        // vector nearestPoint(const label segI, const vector& point) const;
  
        //- Return segment lengths
        tmp<scalarField> segLengths() const;
  
        //- Return midpoint paremeters (arc-length) of spline segments
        tmp<scalarField> midPointParameters() const;
  
        //- Interpolate field from segments to points
        tmp<vectorField> segmentToPointInterpolate
        (
            const vectorField& midPhi
        ) const;  
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
