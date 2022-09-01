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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "setFreePulleyAngularVelocity.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "polyPatchID.H"

#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(setFreePulleyAngularVelocity, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setFreePulleyAngularVelocity,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::setFreePulleyAngularVelocity::setAngularVelocit()
{
    Info << "Writing force and moment history" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beams =
        time_.lookupObject<beamModel>("beamProperties");

    const PtrList<conicalPulleyContactListList> contacts =
        beams.contact().conicalPulleyContacts();

    
    forAll(pulleyIndices_, pI)
    {
        label curPulleyIndex = pulleyIndices_[pI];

        conicalPulley& pulley =
            const_cast<conicalPulley&>(pulleys[curPulleyIndex]);

        label nBeams = beams.nBeams();

        for (label bI=0; bI<nBeams; bI++)
        {
            label nSeg = beams.nCells(bI);

            for (label segI=0; segI<nSeg; segI++)
            {
                
            }
        }
    }

    
    // if (Pstream::master())
    // {
    //     scalar fcn0 =
    //         mag
    //         (
    //             contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //            .normalContactForce()[0]
    //         );

    //     scalar fcn1 =
    //         mag
    //         (
    //             contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //            .normalContactForce()[1]
    //         );

    //     vector T = beams.contact().splines()[beamIndex_].dRdS()[segmentIndex_];
    //     T /= mag(T) + SMALL;

    //     scalarField frictionalContactForce =
    //     (
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //        .transversalFrictionalContactMoment() & T
    //     )/beams.R(beamIndex_);

    //     scalar smallFcn =
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_].smallFcn_;
                        
    //     scalar mf0 =
    //         frictionalContactForce[0]
    //        /(fcn0 + smallFcn);
                            
    //     scalar mf1 =
    //         frictionalContactForce[1]
    //        /(fcn1 + smallFcn);

    //     scalar ttg0 =
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //        .transversalTangentialGap()[0];
        
    //     scalar ttg1 =
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //        .transversalTangentialGap()[1];
        
    //     scalar ttgIncr0 =
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //        .transversalTangentialGapIncrement()[0];
        
    //     scalar ttgIncr1 =
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //        .transversalTangentialGapIncrement()[1];
        
    //     scalar ttgOffset0 =
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //        .transversalTangentialGapOffset()[0];
        
    //     scalar ttgOffset1 =
    //         contacts[beamIndex_][segmentIndex_][pulleyIndex_]
    //        .transversalTangentialGapOffset()[1];

    //     label w = 15;
        
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << mesh.time().value();
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << fcn0;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << fcn1;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << mf0;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << mf1;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << ttg0;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << ttg1;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << ttgIncr0;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << ttgIncr1;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << ttgOffset0;
    //     historyFilePtr_().width(w);
    //     historyFilePtr_() << ttgOffset1 << endl;
                
    //     // historyFilePtr_()
    //     //     << mesh.time().value()
    //     //     << tab << fcn0
    //     //     << tab << fcn1
    //     //     << tab << mf0
    //     //     << tab << mf1
    //     //     << tab << ttg0
    //     //     << tab << ttg1
    //     //     << tab << ttgIncr0
    //     //     << tab << ttgIncr1
    //     //     << tab << ttgOffset0
    //     //     << tab << ttgOffset1
    //     //     << endl;
    // }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::setFreePulleyAngularVelocity::setFreePulleyAngularVelocity
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    pulleyIndices_(dict.lookup("pulleyIndices")),
    changeFrictionTime_(readScalar(dict.lookup("changeFrictionTime")))
{
    Info << "Creating functio object " << name_ << "\n" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    // const fvMesh& mesh =
    //     time_.lookupObject<fvMesh>(regionName_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::setFreePulleyAngularVelocity::start()
{
    // return writeData();
    return true;
}


#if FOAMEXTEND > 40
bool Foam::setFreePulleyAngularVelocity::execute(const bool forceWrite)
#else
bool Foam::setFreePulleyAngularVelocity::execute()
#endif
{
    // const fvMesh& mesh =
    //     time_.lookupObject<fvMesh>(regionName_);

    // const beamModel& beams =
    //     time_.lookupObject<beamModel>("beamProperties");

    // const PtrList<conicalPulley>& pulleys =
    //     beams.conicalPulleys();

    // if (mesh.time().value() >= changeFrictionTime_-SMALL)
    // {
    //     forAll(pulleyIndices_, pI)
    //     {
    //         label curPulleyIndex = pulleyIndices_[pI];

    //         conicalPulley& pulley =
    //             const_cast<conicalPulley&>(pulleys[curPulleyIndex]);

    //         pulley.frictionless() = false;

    //         Info << "Switching frictionless parameter to false for pulley: "
    //              << pI << endl;
    //     }
    // }

    return setAngularVelocity();
}


bool Foam::setFreePulleyAngularVelocity::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
