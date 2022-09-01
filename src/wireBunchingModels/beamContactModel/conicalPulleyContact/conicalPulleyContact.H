/*---------------------------------------------------------------------------* \
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
    Foam::conicalPulleyContact

Description
    This class describes ...

\*---------------------------------------------------------------------------*/

#ifndef conicalPulleyContact_H
#define conicalPulleyContact_H

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

class conicalPulleyContact;

typedef List<conicalPulleyContact> conicalPulleyContactList;
typedef List<conicalPulleyContactList> conicalPulleyContactListList;

    
    
inline Ostream& operator<<(Ostream& os, const conicalPulleyContact& pc);
inline bool operator==
(
    const conicalPulleyContact& a,
    const conicalPulleyContact& b
);

  
/*---------------------------------------------------------------------------*\
                           Class conicalPulleyContact Declaration
\*---------------------------------------------------------------------------*/

class conicalPulleyContact
{
    // Private data

        //- Referenco to list of all contacts
        const PtrList<conicalPulleyContactListList>* conicalPulleyContactsPtr_;
    
        //- Index of first beam in line contact
        mutable label firstBeam_;

        //- First beam segment in line contact
        mutable label firstBeamSegment_;

        //- First beam parameter in line contact
        mutable scalar firstBeamZeta_;

        //- Index of the pulley in line contact
        mutable label pulleyIndex_;

        //- Pulley contact point coordinates
        mutable vectorField pulleyContactCoordinates_;

        //- Pulley contact point coordinates
        mutable vectorField oldPulleyContactCoordinates_;
  
        //- Contact distance
        mutable scalarField delta_;
  
        //- Normal contact force direction
        mutable vectorField normalDir_;
    
        //- Contact angle
        mutable scalarField contactAngle_;
  
        //- Normal contact force
        mutable vectorField normalContactForce_;

        //- Normal contact force
        mutable scalarField normalContactForceOffset_;
    
        //- Normal contact force
        mutable vectorField oldNormalContactForce_;
    
        //- Normal contact force derivative
        mutable tensorField normalContactForceDerivative_;

        //- Normal contact force direction derivative
        mutable tensorField normalContactForceDirectionDerivative_;
    
        //- Normal gap
        mutable scalarField normalGap_;
  
        //- Normal gap
        mutable scalarField oldNormalGap_;
    
        //- Frictional contact moment
        mutable vectorField frictionalContactMoment_;
    
        //- Frictional contact force
        mutable vectorField frictionalContactForce_;
        
        //- Frictional contact moment derivative
        mutable tensorField frictionalContactMomentDerivative_;

        //- Frictional contact moment derivative
        mutable tensorField frictionalContactMomentDerivativeOverFcn_;

        //- Frictional contact force derivative
        mutable tensorField frictionalContactForceDerivative_;

        //- Frictional contact force derivative
        mutable tensorField frictionalContactForceDerivativeOverFcn_;

        //- Torsion component of the rotation angle
        //  after the contact inception
        mutable scalarField torsionTheta_;
        mutable scalarField torsionDTheta_;

        //- Rotation matrix at the old time step
        mutable tensorField oldRM_;


        //- Tangential gap
        mutable scalarField tangentialGap_;
        //- Tangential gap offset
        mutable scalarField tangentialGapOffset_;


        //- Axial tangential gap
        mutable scalarField axialTangentialGap_;

        //- Axial tangential gap
        mutable scalarField axialTangentialGapIncrement_;
    
        //- Axial tangential gap offset
        mutable scalarField axialTangentialGapOffset_;


        //- Axial component of the frictional contact force
        mutable vectorField axialFrictionalContactForce_;

        //- Axial component of the frictional contact force
        mutable tensorField axialFrictionalContactForceDerivative_;


        //- Transversal tangential gap
        mutable scalarField transversalTangentialGap_;

        //- Transversal tangential gap
        mutable scalarField transversalTangentialGapIncrement_;

        //- Transversal tangential gap offset
        mutable scalarField transversalTangentialGapOffset_;

        //- Axial component of the frictional contact force
        mutable vectorField transversalFrictionalContactForce_;

        //- Axial component of the frictional contact force
        mutable tensorField transversalFrictionalContactForceDerivativePerDTheta_;
    
        //- Axial component of the frictional contact force
        mutable tensorField transversalFrictionalContactForceDerivativePerDW_;
        
        //- Axial component of the frictional contact force
        mutable vectorField transversalFrictionalContactMoment_;

        //- Axial component of the frictional contact force
        mutable tensorField transversalFrictionalContactMomentDerivativePerDTheta_;

        //- Axial component of the frictional contact force
        mutable tensorField transversalFrictionalContactMomentDerivativePerDW_;


    
        //- Axial force due to bearing friction
        mutable vectorField axialBearingFrictionForce_;

        //- Normal contact state (0, 1, 2)
        label normalContactState_;
        label oldNormalContactState_;
    
        //- Active friction contact
        mutable List<bool> activeFriction_;

        //- Stick friction contact
        mutable List<bool> stick_;
    
        //- Unloading friction contact
        mutable List<bool> unloading_;
    
        //- Current time index
        mutable label curTimeIndex_;

public:

        static scalar smallFcn_;
    
    // Constructors

        //- Construct null
        conicalPulleyContact();
    
        //- Construct
        conicalPulleyContact
        (
            const PtrList<conicalPulleyContactListList>& conicalPulleyContacts
        );

        //- Construct
        conicalPulleyContact
        (
            const PtrList<conicalPulleyContactListList>& conicalPulleyContacts,
            const label fb,
            const label fbSeg,
            const scalar fbZeta,
            const label pulleyIndex,
            const vectorField& pcc,
	    const label timeIndex
        );

    // Member Functions

        //- Return reference to conical pully contact list
        const PtrList<conicalPulleyContactListList>& conicalPulleyContacts()
        {
            return *conicalPulleyContactsPtr_;
        }
    
        //- Unloading 
        const List<bool>& unloading() const
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
  
        //- Return conical pulley contact point coordinates
        const vectorField& pulleyContactCoordinates() const
        {
            return pulleyContactCoordinates_;
        }
  
        //- Return conical pulley contact point coordinates
        const vectorField& oldPulleyContactCoordinates() const
        {
            return oldPulleyContactCoordinates_;
        }
    
        //- Return point contact distance
        const scalarField& delta() const
        {
            return delta_;
        }
  
        //- Return point contact distance
        scalarField& delta()
        {
            return delta_;
        }

        //- Return current contact force
        const scalarField& normalGap() const
        {
            return normalGap_;
        }
  
        //- Return current contact force
        const vectorField& normalDirection() const
        {
            return normalDir_;
        }
    
        //- Return current contact force
        const vectorField& normalContactForce() const
        {
            return normalContactForce_;
        }
  
        //- Return current contact force
        vectorField& normalContactForce()
        {
            return normalContactForce_;
        }
  
        //- Return current contact force
        const scalarField& normalContactForceOffset() const
        {
            return normalContactForceOffset_;
        }
  
        //- Return current contact force
        scalarField& normalContactForceOffset()
        {
            return normalContactForceOffset_;
        }
    
        //- Return current contact force
        const vectorField& oldNormalContactForce() const
        {
            return oldNormalContactForce_;
        }
    
        //- Return current contact force derivative
        const tensorField& normalContactForceDerivative() const
        {
            return normalContactForceDerivative_;
        }
  
        //- Return current contact force derivative
        tensorField& normalContactForceDerivative()
        {
            return normalContactForceDerivative_;
        }

        //- Return current contact force direction derivative
        const tensorField& normalContactForceDirectionDerivative() const
        {
            return normalContactForceDirectionDerivative_;
        }
  
        //- Return current contact force direction derivative
        tensorField& normalContactForceDirectionDerivative()
        {
            return normalContactForceDirectionDerivative_;
        }
    
        //- Return current contact force
        const vectorField& frictionalContactMoment() const
        {
            return frictionalContactMoment_;
        }

        //- Return current contact force
        const vectorField& frictionalContactForce() const
        {
            return frictionalContactForce_;
        }
    
        //- Return current contact force derivative
        const tensorField& frictionalContactMomentDerivative() const
        {
            return frictionalContactMomentDerivative_;
        }

        //- Return current contact force derivative
        const tensorField& frictionalContactMomentDerivativeOverFcn() const
        {
            return frictionalContactMomentDerivativeOverFcn_;
        }
    
        //- Return current contact force derivative
        const tensorField& frictionalContactForceDerivative() const
        {
            return frictionalContactForceDerivative_;
        }

        //- Return current contact force derivative
        const tensorField& frictionalContactForceDerivativeOverFcn() const
        {
            return frictionalContactForceDerivativeOverFcn_;
        }
    
        //- Return current contact force
        const vectorField& axialBearingFrictionForce() const
        {
            return axialBearingFrictionForce_;
        }

        //- Return current contact force
        vectorField& axialBearingFrictionForce()
        {
            return axialBearingFrictionForce_;
        }

        const scalarField& axialTangentialGap() const
        {
            return axialTangentialGap_;
        }

        const scalarField& axialTangentialGapIncrement() const
        {
            return axialTangentialGapIncrement_;
        }
    
        scalarField& axialTangentialGapIncrement()
        {
            return axialTangentialGapIncrement_;
        }
    
        //- Axial component of frictional contact force
        const vectorField& axialFrictionalContactForce() const
        {
            return axialFrictionalContactForce_;
        }

        //- Derivative of axial component of frictional contact force
        const tensorField& axialFrictionalContactForceDerivative() const
        {
            return axialFrictionalContactForceDerivative_;
        }
    
        //- Transversal component of frictional contact force
        const vectorField& transversalFrictionalContactForce() const
        {
            return transversalFrictionalContactForce_;
        }

        //- Derivative of transversal component of frictional contact force
        const tensorField& transversalFrictionalContactForceDerivativePerDTheta() const
        {
            return transversalFrictionalContactForceDerivativePerDTheta_;
        }

        //- Derivative of transversal component of frictional contact force
        const tensorField& transversalFrictionalContactForceDerivativePerDW() const
        {
            return transversalFrictionalContactForceDerivativePerDW_;
        }

        //- Transversal tangential gap
        const scalarField& transversalTangentialGap() const
        {
            return transversalTangentialGap_;
        }

        //- Transversal tangential gap increment
        const scalarField& transversalTangentialGapIncrement() const
        {
            return transversalTangentialGapIncrement_;
        }

        //- Transversal tangential gap increment
        scalarField& transversalTangentialGapIncrement()
        {
            return transversalTangentialGapIncrement_;
        }
    
        //- Transversal tangential gap offset
        const scalarField& transversalTangentialGapOffset() const
        {
            return transversalTangentialGapOffset_;
        }
    
        //- Transversal component of frictional contact moment
        const vectorField& transversalFrictionalContactMoment() const
        {
            return transversalFrictionalContactMoment_;
        }

        //- Derivative of transversal component of frictional contact moment
        const tensorField& transversalFrictionalContactMomentDerivativePerDTheta() const
        {
            return transversalFrictionalContactMomentDerivativePerDTheta_;
        }

        //- Derivative of transversal component of frictional contact moment
        const tensorField& transversalFrictionalContactMomentDerivativePerDW() const
        {
            return transversalFrictionalContactMomentDerivativePerDW_;
        }

        //- Return active friction
        const List<bool>& activeFriction() const
        {
            return activeFriction_;
        }
    
        void set
        (
	    const label fb,
	    const label fbSeg,
	    const scalar fbZeta,
	    const label pulleyIndex,
	    const vectorField& pulleyContactCoord,
	    const label timeIndex
	);
  
        //- Update normal and tangential contact forces
        void updateNormalForce
        (
            const beamModel& beam,
            const PtrList<HermiteSpline>& splines,
            const normalContactModel& normalModel,
            const frictionContactModel& frictionModel,
	    const scalar contactAngleLimit,
            const bool agumentedLagrangian = false
        );

        //- Update normal and tangential contact forces
        void updateTangentialGap
        (
            const beamModel& beam,
            const PtrList<HermiteSpline>& splines,
            const normalContactModel& normalModel,
            const frictionContactModel& frictionModel,
	    const scalar contactAngleLimit,
            const bool agumentedLagrangian = false
        );
    
        //- Update normal and tangential contact forces
        void updateFrictionForceAndMoment
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
	    const conicalPulleyContact& a,
	    const conicalPulleyContact& b
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

        friend Ostream& operator<<(Ostream& os, const conicalPulleyContact& pc)
        {
             return os << pc.firstBeam() << token::SPACE
                       << pc.firstBeamSegment() << token::SPACE
                       << pc.firstBeamZeta() << token::SPACE
                       << pc.pulleyIndex() << token::SPACE
                       << pc.pulleyContactCoordinates();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
