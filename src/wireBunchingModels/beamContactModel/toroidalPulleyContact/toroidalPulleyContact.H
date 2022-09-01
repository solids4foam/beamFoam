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
    Foam::toroidalPulleyContact

Description
    This class describes ...

\*---------------------------------------------------------------------------*/

#ifndef toroidalPulleyContact_H
#define toroidalPulleyContact_H

#include "bool.H"
#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "HermiteSpline.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
class beamModel;

// Forward declaration of friend functions and operators

class toroidalPulleyContact;
inline Ostream& operator<<(Ostream& os, const toroidalPulleyContact& pc);
inline bool operator==
(
    const toroidalPulleyContact& a,
    const toroidalPulleyContact& b
);

  
/*---------------------------------------------------------------------------*\
                           Class toroidalPulleyContact Declaration
\*---------------------------------------------------------------------------*/

class toroidalPulleyContact
{
    // Private data

        //- Index of first beam in line contact
        mutable label firstBeam_;

        //- First beam segment in line contact
        mutable label firstBeamSegment_;

        //- First beam parameter in line contact
        mutable scalar firstBeamZeta_;

        //- Index of the pulley in line contact
        mutable label pulleyIndex_;

        //- Pulley contact point coordinate
        mutable vector pulleyContactCoordinate_;

        //- Contact distance
        mutable scalar delta_;
  
        //- Normal contact force direction
        mutable vector normalDir_;
    
        //- Contact angle
        mutable scalar contactAngle_;
  
        //- Normal contact force
        mutable vector normalContactForce_;

        //- Normal contact force
        mutable vector oldNormalContactForce_;
    
        //- Normal contact force derivative
        mutable tensor normalContactForceDerivative_;

        //- Normal contact force direction derivative
        mutable tensor normalContactForceDirectionDerivative_;
    
        //- Normal gap
        mutable scalar normalGap_;
  
        //- Normal gap
        mutable scalar oldNormalGap_;
    
        //- Frictional contact moment
        mutable vector frictionalContactMoment_;
    
        //- Frictional contact force
        mutable vector frictionalContactForce_;
    
        //- Frictional contact force due to bearing friction
        mutable vector axialFrictionalContactForce_;
    
        //- Frictional contact moment derivative
        mutable tensor frictionalContactMomentDerivative_;

        //- Frictional contact moment derivative
        mutable tensor frictionalContactMomentDerivativeOverFcn_;
    
        //- Frictional contact force derivative
        mutable tensor frictionalContactForceDerivative_;

        //- Frictional contact force derivative
        mutable tensor frictionalContactForceDerivativeOverFcn_;
        
        //- Torsion component of the rotation angle
        //  after the contact inception
        mutable scalar oldTorsionTheta_;
        mutable scalar torsionDTheta_;

        //- Rotation matrix at the old time step
        mutable tensor oldRM_;

        //- Tangential gap
        mutable scalar tangentialGap_;

        //- Tangential gap offset
        mutable scalar tangentialGapOffset_;

        //- Active friction contact
        mutable bool activeFriction_;

        //- Stick friction contact
        mutable bool stick_;
    
        //- Unloading friction contact
        mutable bool unloading_;
    
        //- Current time index
        mutable label curTimeIndex_;

public:

        static scalar smallFcn_;
    
    // Constructors

        //- Construct null
        toroidalPulleyContact();

        //- Construct
        toroidalPulleyContact
        (
            const label fb,
            const label fbSeg,
            const scalar fbZeta,
            const label pulleyIndex,
            const vector pcc,
	    const label timeIndex
        );

    // Member Functions

        //- Unloading 
        const bool& unloading() const
        {
            return unloading_;
        }
    
        //- Return first beam in line contact
        label firstBeam() const
        {
            return firstBeam_;
        }

        //- Return first beam segment in line contact
        label firstBeamSegment() const
        {
            return firstBeamSegment_;
        }

        //- Return first beam segment parameter in line contact
        scalar firstBeamZeta() const
        {
            return firstBeamZeta_;
        }

        //- Return second beam in point contact
        label pulleyIndex() const
        {
            return pulleyIndex_;
        }
  
        //- Return toroidal pulley contact point coordinates
        const vector& pulleyContactCoordinate() const
        {
            return pulleyContactCoordinate_;
        }
  
        //- Return point contact distance
        const scalar& delta() const
        {
            return delta_;
        }
  
        //- Return point contact distance
        scalar& delta()
        {
            return delta_;
        }

        //- Return current contact force
        const scalar& normalGap() const
        {
            return normalGap_;
        }
  
        //- Return current contact force
        const vector& normalContactForce() const
        {
            return normalContactForce_;
        }
  
        //- Return current contact force
        vector& normalContactForce()
        {
            return normalContactForce_;
        }
  
        //- Return current contact force
        const vector& oldNormalContactForce() const
        {
            return oldNormalContactForce_;
        }
    
        //- Return current contact force derivative
        const tensor& normalContactForceDerivative() const
        {
            return normalContactForceDerivative_;
        }
  
        //- Return current contact force derivative
        tensor& normalContactForceDerivative()
        {
            return normalContactForceDerivative_;
        }

        //- Return current contact force direction derivative
        const tensor& normalContactForceDirectionDerivative() const
        {
            return normalContactForceDirectionDerivative_;
        }
  
        //- Return current contact force direction derivative
        tensor& normalContactForceDirectionDerivative()
        {
            return normalContactForceDirectionDerivative_;
        }
    
        //- Return current contact force
        const vector& frictionalContactMoment() const
        {
            return frictionalContactMoment_;
        }

        //- Return current contact force
        const vector& frictionalContactForce() const
        {
            return frictionalContactForce_;
        }

        //- Return current contact force
        const vector& axialFrictionalContactForce() const
        {
            return axialFrictionalContactForce_;
        }

        //- Return current contact force
        vector& axialFrictionalContactForce()
        {
            return axialFrictionalContactForce_;
        }
    
        //- Return current contact force derivative
        const tensor& frictionalContactMomentDerivative() const
        {
            return frictionalContactMomentDerivative_;
        }

        //- Return current contact force derivative
        const tensor& frictionalContactMomentDerivativeOverFcn() const
        {
            return frictionalContactMomentDerivativeOverFcn_;
        }
    
        //- Return current contact force derivative
        const tensor& frictionalContactForceDerivative() const
        {
            return frictionalContactForceDerivative_;
        }

        //- Return current contact force derivative
        const tensor& frictionalContactForceDerivativeOverFcn() const
        {
            return frictionalContactForceDerivativeOverFcn_;
        }
    
        //- Return active friction
        const bool& activeFriction() const
        {
            return activeFriction_;
        }
    
        void set
        (
	    const label fb,
	    const label fbSeg,
	    const scalar fbZeta,
	    const label pulleyIndex,
	    const vector pulleyContactCoord,
	    const label timeIndex
	);
  
        //- Update normal and tangential contact forces
        void updateForce
        (
            const beamModel& beam,
            const PtrList<HermiteSpline>& splines,
            const normalContactModel& normalModel,
            const frictionContactModel& frictionModel,
	    const scalar contactAngleLimit,
            const bool agumentedLagrangian = false
        );

    // Friend Operators

        friend bool operator==
	(
	    const toroidalPulleyContact& a,
	    const toroidalPulleyContact& b
	)
        {
            return
            (
                (a.firstBeam_ == b.firstBeam_)
             && (a.firstBeamSegment_ == b.firstBeamSegment_)
             && (a.pulleyIndex_ == b.pulleyIndex_)
             // && (a.pulleyContactCoordinates_ == b.pulleyContactCoordinates_)
            );
        }

    // Ostream operator

        friend Ostream& operator<<(Ostream& os, const toroidalPulleyContact& pc)
        {
             return os << pc.firstBeam() << token::SPACE
                       << pc.firstBeamSegment() << token::SPACE
                       << pc.firstBeamZeta() << token::SPACE
                       << pc.pulleyIndex() << token::SPACE
                       << pc.pulleyContactCoordinate();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
