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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "conicalPulleyContactHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "polyPatchID.H"

#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conicalPulleyContactHistory, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        conicalPulleyContactHistory,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::conicalPulleyContactHistory::writeData()
{
    Info << "Writing force and moment history" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);
    
    const beamModel& beams =
        time_.lookupObject<beamModel>("beamProperties");

    const PtrList<conicalPulleyContactListList> contacts =
        beams.contact().conicalPulleyContacts();
    
    if (Pstream::master())
    {
        scalar fcn0 =
            mag
            (
                contacts[beamIndex_][segmentIndex_][pulleyIndex_]
               .normalContactForce()[0]
            );

        scalar fcn1 =
            mag
            (
                contacts[beamIndex_][segmentIndex_][pulleyIndex_]
               .normalContactForce()[1]
            );

        vector T = beams.contact().splines()[beamIndex_].dRdS()[segmentIndex_];
        T /= mag(T) + SMALL;

        scalarField frictionalContactForce =
        (
            contacts[beamIndex_][segmentIndex_][pulleyIndex_]
           .transversalFrictionalContactMoment() & T
        )/beams.R(beamIndex_);

        scalar smallFcn =
            contacts[beamIndex_][segmentIndex_][pulleyIndex_].smallFcn_;
                        
        scalar mf0 =
            frictionalContactForce[0]
           /(fcn0 + smallFcn);
                            
        scalar mf1 =
            frictionalContactForce[1]
           /(fcn1 + smallFcn);

        scalar ttg0 =
            contacts[beamIndex_][segmentIndex_][pulleyIndex_]
           .transversalTangentialGap()[0];
        
        scalar ttg1 =
            contacts[beamIndex_][segmentIndex_][pulleyIndex_]
           .transversalTangentialGap()[1];
        
        scalar ttgIncr0 =
            contacts[beamIndex_][segmentIndex_][pulleyIndex_]
           .transversalTangentialGapIncrement()[0];
        
        scalar ttgIncr1 =
            contacts[beamIndex_][segmentIndex_][pulleyIndex_]
           .transversalTangentialGapIncrement()[1];
        
        scalar ttgOffset0 =
            contacts[beamIndex_][segmentIndex_][pulleyIndex_]
           .transversalTangentialGapOffset()[0];
        
        scalar ttgOffset1 =
            contacts[beamIndex_][segmentIndex_][pulleyIndex_]
           .transversalTangentialGapOffset()[1];

        label w = 15;
        
        historyFilePtr_().width(w);
        historyFilePtr_() << mesh.time().value();
        historyFilePtr_().width(w);
        historyFilePtr_() << fcn0;
        historyFilePtr_().width(w);
        historyFilePtr_() << fcn1;
        historyFilePtr_().width(w);
        historyFilePtr_() << mf0;
        historyFilePtr_().width(w);
        historyFilePtr_() << mf1;
        historyFilePtr_().width(w);
        historyFilePtr_() << ttg0;
        historyFilePtr_().width(w);
        historyFilePtr_() << ttg1;
        historyFilePtr_().width(w);
        historyFilePtr_() << ttgIncr0;
        historyFilePtr_().width(w);
        historyFilePtr_() << ttgIncr1;
        historyFilePtr_().width(w);
        historyFilePtr_() << ttgOffset0;
        historyFilePtr_().width(w);
        historyFilePtr_() << ttgOffset1 << endl;
    }

    if (time_.outputTime())
    {
        OFstream file
        (
            time_.timePath()/"conicalPulleyNormalContactForce.txt"
        );

        label startPatchIndex_ = 0;
        label endPatchIndex_ = 1;
        
        file.precision(12);

        const volVectorField& Fcn =
            beams.contact().conicalPulleyNormalContactForce();
        const vectorField& FcnI = Fcn.internalField();

        volVectorField C
        (
            IOobject
            (
                "Ru",
                time_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("Ru", dimLength, vector::zero)
        );
        C = mesh.C();

        // C += W;
        const vectorField& CI = C.internalField();

        scalar refx = C.boundaryField()[startPatchIndex_][0].x();

        // file << "x" << tab << "y" << tab << "z" << endl; 

        file << C.boundaryField()[startPatchIndex_][0].x() - refx << tab
             << C.boundaryField()[startPatchIndex_][0].y() << tab
             << C.boundaryField()[startPatchIndex_][0].z() << tab
             << Fcn.boundaryField()[startPatchIndex_][0].x() << tab
             << Fcn.boundaryField()[startPatchIndex_][0].y() << tab
             << Fcn.boundaryField()[startPatchIndex_][0].z() << endl;

        forAll(CI, cellI)
        {
            file << CI[cellI].x() - refx << tab
                 << CI[cellI].y() << tab
                 << CI[cellI].z() << tab
                 << FcnI[cellI].x() << tab
                 << FcnI[cellI].y() << tab
                 << FcnI[cellI].z() << endl; 
        }

        file << C.boundaryField()[endPatchIndex_][0].x() - refx << tab
             << C.boundaryField()[endPatchIndex_][0].y() << tab
             << C.boundaryField()[endPatchIndex_][0].z() << tab
             << Fcn.boundaryField()[endPatchIndex_][0].x() << tab
             << Fcn.boundaryField()[endPatchIndex_][0].y() << tab
             << Fcn.boundaryField()[endPatchIndex_][0].z() << endl;

        // Write beam axial force
        {
            OFstream file
            (
                time_.timePath()/"activeNormalContact.txt"
            );

            file.precision(12);

            surfaceVectorField Cf = mesh.Cf();
            
            const volVectorField& Fcn =
                beams.contact().conicalPulleyNormalContactForce();
            const vectorField& FcnI = Fcn.internalField();

            // const surfaceScalarField& Qa =
            //     mesh.lookupObject<surfaceScalarField>("Qa");

            scalar refx = Cf.boundaryField()[startPatchIndex_][0].x();

            scalar curMaxFcn =
                max
                (
                    FcnI[mesh.boundary()[startPatchIndex_].faceCells()[0]].x(),
                    FcnI[mesh.boundary()[startPatchIndex_].faceCells()[0]].y()
                );

            if (curMaxFcn < SMALL)
            {
                curMaxFcn = 0;
            }
            else
            {
                curMaxFcn = 1;
            }
            
            file << Cf.boundaryField()[startPatchIndex_][0].x() - refx << tab
                 << curMaxFcn << endl;

            const labelList& own = mesh.owner();
            const labelList& nei = mesh.neighbour();
            
            forAll(Cf, faceI)
            {
                scalar totalFcn0 =
                    FcnI[own[faceI]].x() +
                    FcnI[nei[faceI]].x();

                scalar totalFcn1 =
                    FcnI[own[faceI]].y() +
                    FcnI[nei[faceI]].y();

                if (totalFcn0 < SMALL)
                {
                    totalFcn0 = 0;
                }
                else
                {
                    totalFcn0 = 1;
                }

                if (totalFcn1 < SMALL)
                {
                    totalFcn1 = 0;
                }
                else
                {
                    totalFcn1 = 1;
                }
                
                file << Cf[faceI].x() - refx << tab
                     << (totalFcn0 + totalFcn1)/2 << endl;
            }

            curMaxFcn =
                max
                (
                    FcnI[mesh.boundary()[endPatchIndex_].faceCells()[0]].x(),
                    FcnI[mesh.boundary()[endPatchIndex_].faceCells()[0]].y()
                );

            if (curMaxFcn < SMALL)
            {
                curMaxFcn = 0;
            }
            else
            {
                curMaxFcn = 1;
            }

            file << Cf.boundaryField()[endPatchIndex_][0].x() - refx << tab
                 << curMaxFcn << endl;
        }        
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conicalPulleyContactHistory::conicalPulleyContactHistory
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
    historyFilePtr_(NULL),
    pulleyIndex_(readInt(dict.lookup("pulleyIndex"))),
    beamIndex_(readInt(dict.lookup("beamIndex"))),
    segmentIndex_(readInt(dict.lookup("segmentIndex")))
{
    Info << "Creating functio object " << name_ << "\n" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Create history file if not already created
    if (historyFilePtr_.empty())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                mesh.time().timeName(mesh.time().startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"history"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"history"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            fileName file("conicalPulleyContactData.txt");

            // Open new file at start up
            historyFilePtr_.reset(new OFstream(historyDir/file));

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                label w = 15;
                
                historyFilePtr_().width(w);
                historyFilePtr_() << "# Time";
                historyFilePtr_().width(w);
                historyFilePtr_() << "fcn0";
                historyFilePtr_().width(w);
                historyFilePtr_() << "fcn1";
                historyFilePtr_().width(w);
                historyFilePtr_() << "mf0";
                historyFilePtr_().width(w);
                historyFilePtr_() << "mf1";
                historyFilePtr_().width(w);
                historyFilePtr_() << "ttg0";
                historyFilePtr_().width(w);
                historyFilePtr_() << "ttg1";
                historyFilePtr_().width(w);
                historyFilePtr_() << "ttgIncr0";
                historyFilePtr_().width(w);
                historyFilePtr_() << "ttgIncr1";
                historyFilePtr_().width(w);
                historyFilePtr_() << "ttgOffset0";
                historyFilePtr_().width(w);
                historyFilePtr_() << "ttgOffset1" << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::conicalPulleyContactHistory::start()
{
    return writeData();
}


#if FOAMEXTEND > 40
bool Foam::conicalPulleyContactHistory::execute(const bool forceWrite)
#else
bool Foam::conicalPulleyContactHistory::execute()
#endif
{
    return writeData();
}


bool Foam::conicalPulleyContactHistory::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
