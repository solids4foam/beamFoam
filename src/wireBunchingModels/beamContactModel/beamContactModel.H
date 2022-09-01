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
    beamContactModel

Description
    Beam bunching model class

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    beamContactModel.C

\*---------------------------------------------------------------------------*/

#ifndef beamContactModel_H
#define beamContactModel_H

#include "IOdictionary.H"
#include "vectorList.H"
#include "autoPtr.H"
#include "cubicSpline.H"
#include "HermiteSpline.H"
#include "labelScalar.H"
#include "labelScalarLabelScalar.H"
#include "volFields.H"

#include "allAngleContact.H"
#include "pointContact.H"
#include "lineContact.H"
#include "conicalPulleyContact.H"
#include "toroidalPulleyContact.H"

#include "normalContactModel.H"
#include "frictionContactModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class beamModel;
  
/*---------------------------------------------------------------------------*\
                           Class beamContactModel Declaration
\*---------------------------------------------------------------------------*/

typedef List<tensor> tensorList;
typedef List<tensorList> tensorListList;

typedef List<lineContact> lineContactList;
typedef List<lineContactList> lineContactListList;

//   typedef List<allAngleContact> allAngleContactList;
//   typedef List<allAngleContactList> allAngleContactListList;

// typedef List<conicalPulleyContact> conicalPulleyContactList;
// typedef List<conicalPulleyContactList> conicalPulleyContactListList;

typedef List<toroidalPulleyContact> toroidalPulleyContactList;
typedef List<toroidalPulleyContactList> toroidalPulleyContactListList;

typedef Pair<scalar> scalarPair;
typedef List<scalarPair> scalarPairList;

class beamContactModel
:
    public IOdictionary
{
    // Private data

        //- All beams in one model
        const beamModel& beam_;

        //- Beams repsented by cubic splines
        mutable PtrList<HermiteSpline> splines_;
    
        //- Normal contact model pointer
        autoPtr<normalContactModel> normalModel_;

        //- Friction contact model
        autoPtr<frictionContactModel> frictionModel_;
	
	//- Normal contact model pointer for point contact
	autoPtr<normalContactModel> pointNormalModel_;
	
	//- Friction contact model for point contact
	autoPtr<frictionContactModel> pointFrictionModel_;

        //- Line contacts
        mutable PtrList<lineContactListList> lineContacts_;

        //- Conical pulley contacts
        mutable PtrList<conicalPulleyContactListList> conicalPulleyContacts_;

        //- Toroidal pulley contacts
        mutable PtrList<toroidalPulleyContactListList> toroidalPulleyContacts_;
    
        //- Conical pulley load
        mutable vectorField conicalPulleyLoads_;

        //- Conical pulley contact gaps
        mutable PtrList<scalarPairList> conicalPulleyContactGaps_;

        //- Active conical pulley contact
        mutable volScalarField activeConicalPulleyContact_;

        //- Conical pulley normal contact force
        mutable volVectorField conicalPulleyNormalContactForce_;
    
        //- Conical pulley normal contact force
        mutable volVectorField conicalPulleyAxialContactForce_;
    
        //- Conical pulley frictional-to-normal contact force ratio
        mutable volVectorField conicalPulleyContactForceRatio_;

        //- Toroidal pulley load
        mutable vectorField toroidalPulleyLoads_;

        //- Toroidal pulley contact gaps
        mutable PtrList<scalarPairList> toroidalPulleyContactGaps_;

        //- Active toroidal pulley contact
        mutable volScalarField activeToroidalPulleyContact_;

        //- Toroidal pulley normal contact force
        mutable volVectorField toroidalPulleyNormalContactForce_;
    
        //- Toroidal pulley frictional-to-normal contact force ratio
        mutable volVectorField toroidalPulleyContactForceRatio_;
    
        //- Nearest distance to neighboring beams for every beam
        mutable PtrList<labelScalarListList> nearestPoints_;

        //- Contact force between beams
        mutable PtrList<vectorListList> contactForces_;

        //- Contact force between beams
        mutable PtrList<tensorListList> contactForceDerivatives_;

        //- Contact angle between current beam and all neighboring beams
        mutable PtrList<scalarListList> contactAngles_;

        //- List of point contacts
        mutable DynamicList<pointContact> pointContacts_;
  
        //- Contact force between beams
        mutable PtrList<vectorListList> frictionalContactForces_;
  
        //- Contact force between beams
        mutable PtrList<vectorListList> frictionalContactMoments_;
  
        //- Gap between beams
        mutable PtrList<scalarListList> contactGaps_;
  
        //- Gap offset
        mutable PtrList<scalarListList> contactOffsets_;

        //- Distance between beams
        mutable PtrList<scalarListList> contactDistances_;
  
        //- Implicit contact
        bool implicit_;

        //- Augmented Lagrangian contact model
        bool augmentedLagrangian_;
  
        //- Implicit contact
        bool sendAllPointsAndTangents_;
    
        //- Contact angle limit ( angles less than this -- line contact)
        scalar lowerContactAngleLimit_;
	
	//- Contact angle upper limit (angles more than this -- point contact)
	scalar upperContactAngleLimit_;
  
        //- Current time index
        label curTimeIndex_;

        //- Total splines update time
        scalar totalSplinesUpdateTime_;
    
    // Private Member Functions

        //- Initialize line contacts
        void initializeLineContacts();

        //- Initialize line contact
        void initializeLineContact
		(
            const label bI,
			const label segI,
			const label nbI
		);

        //- Update line contacts
        void updateLineContacts();

        //- Update line contact
        void updateLineContact
		(
            const label bI,
			const label segI,
			const label nbI
		);

        //- Update line contact forces
        void updateLineContactForces();

        //- Update normal gap offset for penalty contact law
        void updateNormalGapOffset();

        //- Update conical pulley contacts
        void updateConicalPulleyContacts();

        //- Update conical pulley contact forces
        void updateConicalPulleyContactForces();

        //- Update toroidal pulley contacts
        void updateToroidalPulleyContacts();

        //- Update toroidal pulley contact forces
        void updateToroidalPulleyContactForces();

        //- Initialize point contacts
        void initializePointContacts();

        //- Add new point contacts
        void addNewPointContacts();

        //- Update existing point contacts
        void updatePointContacts();

        //- Update point contact forces
        void updatePointContactForces();

        //- Disallow default bitwise copy construct
        beamContactModel(const beamContactModel&);

        //- Disallow default bitwise assignment
        void operator=(const beamContactModel&);

protected:

    // Protected member functions

public:

    //- Runtime type information
    TypeName("beamContactModel");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from components
        beamContactModel
        (
            const beamModel& beam
        );

    // Selectors


    // Destructor

        virtual ~beamContactModel();

    // Member Functions

        // Access

            //- Return number of contact correctors
            bool implicit() const
            {
                return implicit_;
            }

            //- Check if contact is active
            bool active() const;

            //- Return beams represented by cubic splines
            const PtrList<HermiteSpline>& splines() const
            {
                return splines_;
            }

            //- Normal contact model pointer
            const normalContactModel& normalModel()
            {
                return normalModel_();
            }
    
            //- Return updated nearest points
            const PtrList<conicalPulleyContactListList>&
            conicalPulleyContacts() const
            {
                return conicalPulleyContacts_;
            }

            //- Return updated nearest points
            PtrList<conicalPulleyContactListList>&
            conicalPulleyContacts()
            {
                return conicalPulleyContacts_;
            }
    
            //- Return conical pulley loads
            const vectorField& conicalPulleyLoads() const
            {
                return conicalPulleyLoads_;
            }
    
            //- Return conical pulley contact gaps
            const PtrList<scalarPairList>&
            conicalPulleyContactGaps() const
            {
                return conicalPulleyContactGaps_;
            }
    
            //- Return conical pulley contat indicator
            const volScalarField& activeConicalPulleyContact() const
            {
                return activeConicalPulleyContact_;
            }

            //- Return conical pulley contat indicator
            const volVectorField& conicalPulleyNormalContactForce() const
            {
                return conicalPulleyNormalContactForce_;
            }

            //- Return updated nearest points
            const PtrList<toroidalPulleyContactListList>&
            toroidalPulleyContacts() const
            {
                return toroidalPulleyContacts_;
            }

            //- Return toroidal pulley loads
            const vectorField& toroidalPulleyLoads() const
            {
                return toroidalPulleyLoads_;
            }
    
            //- Return toroidal pulley contact gaps
            const PtrList<scalarPairList>&
            toroidalPulleyContactGaps() const
            {
                return toroidalPulleyContactGaps_;
            }
    
            //- Return toroidal pulley contat indicator
            const volScalarField& activeToroidalPulleyContact() const
            {
                return activeToroidalPulleyContact_;
            }

    
            //- Return updated nearest points
            const PtrList<lineContactListList>& lineContacts() const
            {
                return lineContacts_;
            }
  
            //- Return updated nearest points
            const PtrList<labelScalarListList>& nearestPoints() const
            {
                return nearestPoints_;
            }

            //- Return updated nearest points
            const PtrList<scalarListList>& contactAngles() const
            {
                return contactAngles_;
            }
    
            //- Return updated contact forces
            const PtrList<vectorListList>& contactForces() const
            {
                return contactForces_;
            }
  
            //- Return updated contact force derivatives
            const PtrList<tensorListList>& contactForceDerivatives() const
            {
                return contactForceDerivatives_;
            }

            const DynamicList<pointContact>& pointContacts() const
            {
                return pointContacts_;
            }

            //- Return updated contact forces
            const PtrList<scalarListList>& contactGaps() const
            {
                return contactGaps_;
            }
  
            //- Return updated contact forces
            const PtrList<scalarListList>& contactDistances() const
            {
                return contactDistances_;
            }
  
        // Edit

            //- Update contact (calc. contact forces and their derivaatives)
            scalar update();

            //- Update at the end of time step calculation
            bool finalUpdate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
