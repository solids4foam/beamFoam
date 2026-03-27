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
    Foam::allAngleContact

Description
    This class describes a (contact angel dependent) combination of ponit contact and line-to-line contact
    * It is invoked when the contact angle between two beams is between two defined contact angles alpha1 and alpha2. For, angles less than alpha1, line contact is fully adopted and for angles more than alpha2, only point contact is used .
    

\*---------------------------------------------------------------------------*/

#ifndef allAngleContact_H
#define allAngleContact_H

#include "bool.H"
#include "label.H"
#include "scalar.H"
#include "vector.H"
#include "tensor.H"
#include "HermiteSpline.H"
#include "normalContactModel.H"
#include "frictionContactModel.H"

#include "pointContact.H"
#include "lineContact.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  
class beamModel;

// Forward declaration of friend functions and operators

class allAngleContact;
inline Ostream& operator<<(Ostream& os, const allAngleContact& pc);
inline bool operator==(const allAngleContact& a, const allAngleContact& b);

  
/*---------------------------------------------------------------------------*\
                           Class lineContact Declaration
\*---------------------------------------------------------------------------*/

class allAngleContact
{
    // Private data

        //- Index of first beam in point contact
        mutable label firstBeam_;

        //- First beam segment in point contact
        mutable label firstBeamSegment_;

        //- First beam segment in point contact
        mutable label firstBeamOldSegment_;
  
        //- First beam parameter in point contact
        mutable scalar firstBeamZeta_;

        //- First beam parameter in point contact
        mutable scalar firstBeamOldZeta_;

        //- First beam parameter in point contact
        mutable vector contactDirection_;

        //- First beam parameter in point contact
        mutable vector oldContactDirection_;
    
        //- Index of second beam in point contact
        mutable label secondBeam_;

        //- Second beam segment in point contact
        mutable label secondBeamSegment_;

        //- Second beam lower segment
        mutable label secondBeamLowerSegment_;

        //- Second beam lower segment
        mutable label secondBeamUpperSegment_;

        //- Second beam weight
        mutable scalar secondBeamWeight_;
    
        //- Second beam segment in point contact
        mutable label secondBeamOldSegment_;

        //- Second beam parameter in point contact
        mutable scalar secondBeamZeta_;
    
        // //- Second beam parameter in point contact
        // mutable scalar initialSecondBeamZeta_;
  
        //- Second beam parameter in point contact
        mutable scalar secondBeamOldZeta_;

        //- Contact distance
        mutable scalar delta_;
  
        //- Contact angle
        mutable scalar contactAngle_;
	
        //- Normal contact force
        mutable vector normalContactForce_;

        //- Normal contact force
        mutable vector prevNormalContactForce_;
    
        //- Normal contact force offset
        mutable scalar normalContactForceOffset_;

        //- Normal contact force derivative
        mutable tensor normalContactForceDerivative_;
  
        //- Tangential contact force
        mutable vector axialContactForce_;
  
        //- Tangetial contact force derivative
        mutable tensor axialContactForceDerivative_;
  
        //- Tangential contact force
        mutable vector circumferentialContactForce_;
  
        //- Tangential contact force couple
        mutable vector circumferentialContactMoment_;
  
        //- Tangetial contact force derivative
        mutable tensor circumferentialContactForceDerivative_;

        //- Tangetial contact force derivative
        mutable tensor circumferentialContactForceThetaDerivative_;
  
        //- Tangetial contact force derivative
        mutable tensor circumferentialContactMomentDerivative_;
  
        //- Normal gap
        mutable scalar normalGap_;

        //- Normal gap offset
        mutable scalar normalGapOffset_;

        //- Total tangential gap
        mutable scalar axialGap_;

        //- Total tangential gap
        mutable scalar oldAxialGap_;

        //- Total tangential gap
        mutable scalar circumferentialGap_;
  
        //- Total tangential gap
        mutable scalar oldCircumferentialGap_;

        //- Slip gap
        mutable scalar axialSlipGap_;

        //- Slip gap
        mutable scalar oldAxialSlipGap_;

        //- Slip gap
        mutable scalar circumferentialSlipGap_;

        //- Slip gap
        mutable scalar oldCircumferentialSlipGap_;

        //- Current time index
        label curTimeIndex_;
  
public:

    // Constructors

        //- Construct null
        allAngleContact();

        //- Construct
        allAngleContact
        (
            const label fb,
            const label fbSeg,
            const scalar fbZeta,
            const label sb,
            const label sbSeg,
            const scalar sbZeta,
	    const label timeIndex
        );

    // Member Functions
    
	//- Return line contact normal force
	vector lineContactNormalForce() const
	{
	    // vector i(1,0,0);
	    return lineContact().normalContactForce();
	}
	
	//- Return point contact normal force
	vector pointContactNormalForce() const
	{
	    // vector i(1,0,0);
	    return pointContact().normalContactForce();
	}

        //- Return first beam in point contact
        label firstBeam() const
        {
            return firstBeam_;
        }

        //- Return first beam segment in point contact
        label firstBeamSegment() const
        {
            return firstBeamSegment_;
        }

        //- Return first beam segment parameter in point contact
        scalar firstBeamZeta() const
        {
            return firstBeamZeta_;
        }

        //- Return second beam in point contact
        label secondBeam() const
        {
            return secondBeam_;
        }
  
        //- Return second beam segment in point contact
        label secondBeamSegment() const
        {
            return secondBeamSegment_;
        }

        //- Return second beam lower segment
        label secondBeamLowerSegment() const
        {
            return secondBeamLowerSegment_;
        }

        //- Return second beam upper segment
        label secondBeamUpperSegment() const
        {
            return secondBeamUpperSegment_;
        }
    
        //- Return second beam weight
        scalar secondBeamWeight() const
        {
            return secondBeamWeight_;
        }
    
        //- Return second beam segment parameter in point contact
        scalar secondBeamZeta() const
        {
            return secondBeamZeta_;
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

        //- Return normal gap offset in penalty law
        const scalar& normalGapOffset() const
        {
            return normalGapOffset_;
        }

        //- Return normal gap offset in penalty law
        scalar& normalGapOffset()
        {
            return normalGapOffset_;
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
  
        //- Return current contact force offset
        const scalar& normalContactForceOffset() const
        {
            return normalContactForceOffset_;
        }
  
        //- Return current contact force offset
        scalar& normalContactForceOffset()
        {
            return normalContactForceOffset_;
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

        //- Return current tangential contact force
        const vector& axialContactForce() const
        {
            return axialContactForce_;
        }
  
        //- Return current tangential contact force derivative
        const tensor& axialContactForceDerivative() const
        {
            return axialContactForceDerivative_;
        }
  
        //- Return current tangential contact force
        const vector& circumferentialContactForce() const
        {
            return circumferentialContactForce_;
        }
  
        //- Return current tangential contact moment
        const vector& circumferentialContactMoment() const
        {
            return circumferentialContactMoment_;
        }
  
        //- Return current tangential contact force derivative
        const tensor& circumferentialContactForceDerivative() const
        {
            return circumferentialContactForceDerivative_;
        }
  
        //- Return current tangential contact force derivative
        const tensor& circumferentialContactForceThetaDerivative() const
        {
            return circumferentialContactForceThetaDerivative_;
        }
  
        //- Return current tangential contact force derivative
        const tensor& circumferentialContactMomentDerivative() const
        {
            return circumferentialContactMomentDerivative_;
        }
  
        void set
        (
	    const label fb,
	    const label fbSeg,
	    const scalar fbZeta,
	    const label sb,
	    const label sbSeg,
	    const scalar sbZeta,
	    const label timeIndex
	);
  
        //- Update normal and tangential contact forces
        void updateForce
        (
            const beamModel& beam,
            const PtrList<HermiteSpline>& splines,
            const normalContactModel& normalModel,
            const frictionContactModel& frictionModel,
	    const scalar lowerContactAngleLimit,
	    const scalar upperContactAngleLimit,
            const bool agumentedLagrangian = false
        );
	
	//- Return all angle point contact weight value for 
	//  combining point and line contact
	scalar allAnglePointContactWeight
	(
	    const beamModel& beam,
            const PtrList<HermiteSpline>& splines,
	    const scalar lowerContactAngleLimit,
	    const scalar upperContactAngleLimit	    
	);

	//- Return all angle line contact weight value for 
	//  combining point and line contact
	scalar allAngleLineContactWeight
	(
	    const beamModel& beam,
            const PtrList<HermiteSpline>& splines,
	    const scalar lowerContactAngleLimit,
	    const scalar upperContactAngleLimit	    
	);

        //- Calculate and return current normal gap
        scalar calcNormalGap
        (
            const beamModel& beam,
            const PtrList<HermiteSpline>& splines
        ) const;
    
    // Friend Operators

        friend bool operator==(const allAngleContact& a, const allAngleContact& b)
        {
            return
            (
                (a.firstBeam_ == b.firstBeam_)
             && (a.firstBeamSegment_ == b.firstBeamSegment_)
             && (a.secondBeam_ == b.secondBeam_)
             && (a.secondBeamSegment_ == b.secondBeamSegment_)
            );
        }

    // Ostream operator

        friend Ostream& operator<<(Ostream& os, const allAngleContact& pc)
        {
             return os << pc.firstBeam() << token::SPACE
                       << pc.firstBeamSegment() << token::SPACE
                       << pc.firstBeamZeta() << token::SPACE
                       << pc.secondBeam() << token::SPACE
                       << pc.secondBeamSegment() << token::SPACE
                       << pc.secondBeamZeta();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
