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
    beamsBunchingModel

Description
    Beam bunching model class

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    beamsBunchingModel.C

\*---------------------------------------------------------------------------*/

#ifndef beamsBunchingModel_H
#define beamsBunchingModel_H

#include "beamModel.H"
#include "IOdictionary.H"
#include "Tuple2.H"
#include "vectorList.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class beamsBunchingModel Declaration
\*---------------------------------------------------------------------------*/

typedef Tuple2<label, scalar> labelScalar;
typedef List<labelScalar> labelScalarList;
typedef List<labelScalarList> labelScalarListList;

class beamsBunchingModel
:
    public IOdictionary
{
    // Private data

        //- Reference to the master time object
        Time& runTime_;

        //- All beams
        mutable PtrList<beamModel> beams_;

        //- Nearest distance to neighboring beams for every beam
        mutable PtrList<labelScalarListList> nearestPoints_;
        // mutable PtrList<vectorListList> nearestPoints_;
 
        //- Contact force between beams
        mutable PtrList<vectorListList> contactForces_;

        //- Contact force between beams
        mutable PtrList<vectorListList> frictionalContactForces_;
  
        //- Contact force between beams
        mutable PtrList<vectorListList> frictionalContactMoments_;
  
        //- Gap between beams
        mutable PtrList<scalarListList> contactGap_;
  
        //- Gap offset
        mutable PtrList<scalarListList> contactOffset_;

        //- Friction coefficient
        scalar frictionCoeff_;
  
        //- Contact force relaxation factor
        scalar forceRelaxFac_;

        //- Gap relaxation factor
        scalar gapRelaxFac_;
  
        //- Convergence tolerance
        scalar convergenceTol_;

        //- Number of contact iterations
        label nCorr_;
  
        //- Residual file
        autoPtr<OFstream> residualFilePtr_;

        //- Check if contact is active
        bool active_;
  
    // Private Member Functions

        //- Return contact force as a function of gap 
        scalar contactForce(const scalar g) const;
  
        //- Update contact (calc. gaps and contact forces)
        scalar updateContact();

        //- Update contact forces
        void setContactForces();

        //- Disallow default bitwise copy construct
        beamsBunchingModel(const beamsBunchingModel&);

        //- Disallow default bitwise assignment
        void operator=(const beamsBunchingModel&);
 
protected:

    // Protected member functions

public:

    //- Runtime type information
    TypeName("beamsBunchingModel");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from components
        beamsBunchingModel
        (
            Time& runTime
        );

    // Selectors


    // Destructor

        virtual ~beamsBunchingModel();

    // Member Functions

        // Access

            //- Check if contact is active
            bool active() const
            {
                return active_;
            }
  
            //- Return time
            const Time& runTime() const
            {
                return runTime_;
            }

            //- Return time
            Time& runTime()
            {
                return runTime_;
            }
  
        // Edit

            //- Solve bunching
            void evolve();
  
            //- Update total fields for all beams
            void updateTotalFields();
  
            //- Write fields
            void writeFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
