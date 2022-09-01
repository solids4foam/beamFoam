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
    Foam::toroidalPulley

Description

SourceFiles
    toroidalPulley.C

\*---------------------------------------------------------------------------*/

#ifndef toroidalPulley_H
#define toroidalPulley_H

#include "vector.H"
#include "tensor.H"
#include "dictionary.H"
#include "transformField.H"
#include "triSurface.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
// class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class toroidalPulley;
// Istream& operator>>(Istream&, toroidalPulley&);
inline Ostream& operator<<(Ostream&, const toroidalPulley&);


/*---------------------------------------------------------------------------*\
                          Class toroidalPulley Declaration
\*---------------------------------------------------------------------------*/

class toroidalPulley
{
    // Private data

        //- Origin of the pulley
        vector origin_;

        //- Shift of the pulley origin
        vector originShift_;
    
        //- Temporal displacement of the pulley origin
        vector temporalOriginDisplacement_;
    
        //- Axis (longitudinal) of the pulley
        vector axis_;
    
        //- Polar axis of the pulley
        vector polarAxis_;

        //- Minimal radius of pulley
        scalar rMin_;
    
        //- Radius of toroidal surface
        scalar rToroidalSurf_;

        //- Width of pulley
        scalar width_;

        //- Rotation direction
        label rotDir_;

        //- Mean bearing diameter
        scalar meanBearingDiameter_;

        //- Bearing friction coefficient
        scalar bearingFrictionCoeff_;

        //- Origin displacement time series
        interpolationTable<vector> originDispSeries_;
            
        //- Stl file name
        fileName stlFileName_;

        //- Tri surface
        triSurface stlModel_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        toroidalPulley(const toroidalPulley&);

        //- Disallow default bitwise assignment
        void operator=(const toroidalPulley&);

public:

    // Static data members

        // //- Static data staticData
        // static const dataType staticData;

    // Constructors

        //- Construct null
        toroidalPulley();

        //- Construct from components
        toroidalPulley
        (
            const vector& origin,
            const vector& originShift,
            const vector& axis,
            const vector& polarAxis,
            const scalar& rMin,
            const scalar& rToroidalSurf,
            const scalar& width,
            const label& rotDir
        );

        //- Construct from dictionary
        toroidalPulley(const dictionary& dict);

        // //- Construct from Istream
        // toroidalPulley(Istream&);

        // //- Construct as copy
        // toroidalPulley(const toroidalPulley&);

    // Selectors

        // //- Select null constructed
        // static autoPtr<toroidalPulley> New();

    //- Destructor
    ~toroidalPulley();

    // Member Functions

        // Access

            //- Return origin
            const vector origin() const
            {
                return origin_ + originShift_ + temporalOriginDisplacement_;
            }
    
            //- Return temporal origin displacement
            const vector& temporalOriginDisplacement() const
            {
                return temporalOriginDisplacement_;
            }
    
            //- Return axis
            const vector& axis() const
            {
                return axis_;
            }
    
            //- Return polar axis
            const vector& polarAxis() const
            {
                return polarAxis_;
            }
    
            //- Minimal radius
            scalar minRadius() const
            {
                return rMin_;
            }

            //- Maximal radius
            scalar maxRadius() const
            {
                if (rToroidalSurf_/rMin_ > 1000)
                {
                    return rMin_;
                }
                else
                {
                    return rMin_ + rToroidalSurf_
                      - ::sqrt(sqr(rToroidalSurf_) - sqr(width()/2));
                }
            }
    
            //- Width
            scalar width() const
            {
                return width_;
            }
    
            //- Rotation direction
            const label& rotDir() const
            {
                return rotDir_;
            }

            //- Mean bearing diameter
            const scalar& meanBearingDiameter() const
            {
                return meanBearingDiameter_;
            }

            //- Bearing friction coefficient
            const scalar& bearingFrictionCoeff() const
            {
                return bearingFrictionCoeff_;
            }
    
            //- Nearest point
            vector nearestPoint(const vector& p) const;

            //- Cartesion coordinates for given cylindrical coordinates
            vector position(const vector& cylCoord) const;

            //- Axial tangent vector in cartesian coordinates
            vector axialTangent(const vector& cylCoord) const;
    
            //- dr/dz
            scalar DrDz(const vector& cylCoord) const;
    
            //- Return STL model
            const triSurface& stlModel() const
            {
                return stlModel_;
            }
   
            //- Return STL model
            triSurface& stlModel()
            {
                return stlModel_;
            }
    
        // Check

            //- Check if pulley is moving
            bool moving() const
            {
                return (originDispSeries_.size()>0);
            }
    
        // Edit

            //- Move pulley origin
            void move(const scalar time);
    
        // Write


    // Member Operators

        // void operator=(const toroidalPulley&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        // friend Istream& operator>>(Istream&, toroidalPulley&);
        friend Ostream& operator<<(Ostream& os, const toroidalPulley& tp)
        {
            vector rotAngle(0, 0, 0);

            vector k(0, 0, 1);

            if ( mag(tp.axis() & k) < (1.0-SMALL) )
            {
                tensor R = rotationTensor(k, tp.axis());

                vector rotAxis
                (
                    R.zy() - R.yz(),
                    R.xz() - R.zx(),
                    R.yx() - R.xy()
                );

                scalar magRotAngle = asin(mag(rotAxis)/2)*180/M_PI;

                rotAxis /= mag(rotAxis);

                rotAngle = magRotAngle*rotAxis;
                
                Info << tp.axis() << ", " << rotAngle << endl;
            }

            os << "translate(["
               << tp.origin().x() << ", "
               << tp.origin().y() << ", "
               << tp.origin().z() << "])" << endl;

            os << "rotate(["
               << rotAngle.x() << ", "
               << rotAngle.y() << ", "
               << rotAngle.z() << "])" << endl;

            // os << "        translate(["
            //    << 0 << ", "
            //    << 0 << ", "
            //    << -tp.width()/2 << "])" << endl;

            if (tp.rToroidalSurf_/tp.rMin_ > 1000)
            {
                os << "cylinder("
                   << tp.width() << ", "
                   << tp.minRadius() << ", "
                   << tp.minRadius() << ", "
                   << "true);\n" << endl;
            }
            else
            {
                os << "difference() {";

                os << "cylinder("
                   << tp.width() << ", "
                   << tp.maxRadius() << ", "
                   << tp.maxRadius() << ", "
                   << "true);\n";

                os << "rotate_extrude(convexity=10) ";
                
                os << "translate(["
                   << tp.minRadius() + tp.rToroidalSurf_ << ", "
                   << 0 << ", "
                   << 0 << "]) ";

                os << "circle("
                   << tp.rToroidalSurf_
                   << ");}\n" << endl;
            }
            
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
