/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    plasticityStressResultantReturn

Description
    Method to calcluate the plastic strain increment, returning the
    stress to the yield surface

SourceFiles
    plasticityStressResultantReturn.C
    newPlasticityStressResultantReturn.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef plasticityStressResultantReturn_H
#define plasticityStressResultantReturn_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    
    class beamModel;

/*---------------------------------------------------------------------------*\
                         Class plasticityStressResultantReturn Declaration
\*---------------------------------------------------------------------------*/

class plasticityStressResultantReturn
{
    // Private data

        //- Name
        const word name_;

        //- Plasticity properties dictionary
        dictionary plasticityProperties_;

        //- Current time index
        label curTimeIndex_;

        //- Number of elastic predictor iterations
        label nElPredCorr_;

        //- Elastic prediction tolerance
        scalar elPredTol_;
    
    // Private Member Functions

        //- Disallow copy construct
        plasticityStressResultantReturn(const plasticityStressResultantReturn&);

        //- Disallow default bitwise assignment
        void operator=(const plasticityStressResultantReturn&);


  //protected:

public:

    //- Runtime type information
    TypeName("plasticityStressResultantReturn");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            plasticityStressResultantReturn,
            dictionary,
            (
                const word name,
                beamModel& beamModel
            ),
            (name, beamModel)
        );


    // Selectors

        //- Return a reference to the selected rheology model
        static autoPtr<plasticityStressResultantReturn> New
        (
            const word& name,
            beamModel& beamModel
        );


    // Constructors

        //- Construct from dictionary
        plasticityStressResultantReturn
        (
            const word& name,
            beamModel& beamModel
        );

    // Destructor

        virtual ~plasticityStressResultantReturn()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return current time index
        label& timeIndex()
        {
            return curTimeIndex_;
        }

        //- Return number of elastic prediction correctors
        const label& nElPredCorr()
        {
            return nElPredCorr_;
        }

        //- Return elastic prediction tolerance
        const scalar& elPredTol()
        {
            return elPredTol_;
        }
    
        //- Return flow properties dictionary
        const dictionary& plasticityProperties() const
        {
            return plasticityProperties_;
        }
    
        //- If plasticity is active
        virtual bool plasticityActive() const
        {
            return true;
        }

        //- Increment of plastic strain
        virtual const surfaceVectorField& DGammaP() const = 0;
        virtual const surfaceVectorField& DKP() const = 0;

        //- Elasto-plastic moduli
        virtual const surfaceTensorField& DQDGamma() const = 0;
        virtual const surfaceTensorField& DMDK() const = 0;
        virtual const surfaceTensorField& DMDGamma() const = 0;
        virtual const surfaceTensorField& DQDK() const = 0;

        //- Update the yield stress
        //- called at end of time-step
        virtual void updateYieldStress() = 0;

        //- Correct model
        virtual void correct() = 0;

        //- Write plasticity informations
        virtual void writeInfo() const {}
    
        //- Write model dictionary
        virtual void writeDict(Ostream& os) const {}
    
        //- Write plasticity informations
        virtual void writeFields() const {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
