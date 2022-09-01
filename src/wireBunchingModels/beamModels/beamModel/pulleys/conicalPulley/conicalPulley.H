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
    Foam::conicalPulley

Description

SourceFiles
    conicalPulleyI.H
    conicalPulley.C
    conicalPulleyIO.C

\*---------------------------------------------------------------------------*/

#ifndef conicalPulley_H
#define conicalPulley_H

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
class conicalPulley;

// Istream& operator>>(Istream&, conicalPulley&);
inline Ostream& operator<<(Ostream&, const conicalPulley&);


/*---------------------------------------------------------------------------*\
                          Class conicalPulley Declaration
\*---------------------------------------------------------------------------*/

class conicalPulley
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
    
        //- Angle of pulley (deg)
        scalar angle_;

        //- Width of pulley
        scalar width_;

        //- Rotation direction
        label rotDir_;

        //- Mean bearing diameter
        scalar meanBearingDiameter_;

        //- Bearing friction coefficient
        scalar bearingFrictionCoeff_;

        //- Is it pulley frictionless
        bool frictionless_;
    
        //- Is it pulley frictionless in case when full contact exist
        //  (wire is in contact with both flanges)
        bool frictionlessFullContact_;
    
        //- Apex of the negative side cone
        vector negConeApex_;

        //- Apex of the positive side cone
        vector posConeApex_;

        //- Origin displacement time series
        interpolationTable<vector> originDispSeries_;
            
        //- Stl file name
        fileName stlFileName_;

        //- Tri surface
        triSurface stlModel_;

        //- Circumferential speed at root (m/s) 
        scalar circumferentialSpeed_;
    
        //- Angular velocity (rad/s) 
        mutable scalar angularVelocity_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        conicalPulley(const conicalPulley&);

        //- Disallow default bitwise assignment
        void operator=(const conicalPulley&);

public:

    // Static data members

        // //- Static data staticData
        // static const dataType staticData;

    // Constructors

        //- Construct null
        conicalPulley();

        //- Construct from components
        conicalPulley
        (
            const vector& origin,
            const vector& originShift,
            const vector& axis,
            const vector& polarAxis,
            const scalar& rMin,
            const scalar& angle,
            const scalar& width,
            const label& rotDir
        );

        //- Construct from dictionary
        conicalPulley
        (
            const dictionary& dict
        );

        // //- Construct from Istream
        // conicalPulley(Istream&);

        // //- Construct as copy
        // conicalPulley(const conicalPulley&);

    // Selectors

        // //- Select null constructed
        // static autoPtr<conicalPulley> New();

    //- Destructor
    ~conicalPulley();

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
                return rMin_ + (width_/2)*tan(angle_*M_PI/180);
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

            //- Frictionless pulley
            const bool& frictionless() const
            {
                return frictionless_;
            }

            //- Frictionless pulley
            const bool& frictionlessFullContact() const
            {
                return frictionlessFullContact_;
            }
    
            //- Frictionless pulley
            bool& frictionless()
            {
                return frictionless_;
            }
    
            //- Apex of the negative side cone
            const vector& negConeApex() const
            {
                return negConeApex_;
            }
    
            //- Apex of the positive side cone
            const vector& posConeApex() const
            {
                return posConeApex_;
            }
                
            //- Nearest point at negative cone
            vector nearestPointNeg(const vector& p) const;

            //- Nearest point at positive cone
            vector nearestPointPos(const vector& p) const;

            //- Cartesion coordinates for given cylindrical coordinates
            vector position(const vector& cylCoord) const;

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

            //- Radius of full contact point
            scalar fullContactPointRadius(const scalar wireR) const
            {
                scalar halfAngle = angle_*M_PI/(2*180);
                scalar dr = wireR*cos(halfAngle)/tan(halfAngle);
                scalar contactPointRadius = rMin_ + dr;

                // Info << rMin_ << ", " << dr << endl;
                
                return contactPointRadius;
            }
    
            //- Angular velocity
            const scalar& angularVelocity() const
            {
                if (mag(circumferentialSpeed_) > SMALL)
                {
                    angularVelocity_ =
                        circumferentialSpeed_
                       /minRadius();
                }
                
                return angularVelocity_;
            }
    
            //- Angular velocity
            const scalar& angularVelocity(const scalar wireR) const
            {
                if (mag(circumferentialSpeed_) > SMALL)
                {
                    scalar halfAngle = angle_*M_PI/(2*180);
                    scalar dr = wireR*cos(halfAngle)/tan(halfAngle);
                    
                    scalar contactPointRadius =
                        rMin_ + dr;
                    
                    angularVelocity_ =
                        circumferentialSpeed_
                       /contactPointRadius;
                }
                
                return angularVelocity_;
            }
    
            //- Angular velocity
            void setAngularVelocity(scalar omega)
            {
                if (mag(circumferentialSpeed_) > SMALL)
                {
                    angularVelocity_ =
                        circumferentialSpeed_
                       /minRadius();
                }
                else
                {
                    angularVelocity_ = omega;
                }
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

        // void operator=(const conicalPulley&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        // friend Istream& operator>>(Istream&, conicalPulley&);
        friend Ostream& operator<<(Ostream& os, const conicalPulley& cp)
        {
            vector rotAngle(0, 0, 0);

            vector k(0, 0, 1);

            if ( mag(cp.axis() & k) < (1.0-SMALL) )
            {
                tensor R = rotationTensor(k, cp.axis());

                vector rotAxis
                (
                    R.zy() - R.yz(),
                    R.xz() - R.zx(),
                    R.yx() - R.xy()
                );

                scalar magRotAngle = asin(mag(rotAxis)/2)*180/M_PI;

                rotAxis /= mag(rotAxis);

                rotAngle = magRotAngle*rotAxis;
                
                Info << cp.axis() << ", " << rotAngle << endl;
            }

            os << "translate(["
               << cp.origin().x() << ", "
               << cp.origin().y() << ", "
               << cp.origin().z() << "])" << endl;
            
            os << "    rotate(["
               << rotAngle.x() << ", "
               << rotAngle.y() << ", "
               << rotAngle.z() << "])" << endl;

            os << "        translate(["
               << 0 << ", "
               << 0 << ", "
               << -cp.width()/2 << "])" << endl;
            
            os << "            cylinder("
               << cp.width()/2 << ", "
               << cp.maxRadius() << ", "
               << cp.minRadius() << ", "
               << "false);\n" << endl;

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "conicalPulleyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
