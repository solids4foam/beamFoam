/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "writeBeamData.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "boundBox.H"

#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeBeamData, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeBeamData,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::writeBeamData::writeBeamData
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion)
{
    Info << "Creating writeBeamData function object" << endl;
    
    if (Pstream::parRun())
    {
        FatalErrorIn("writeBeamData::writeBeamData(...)")
            << "writeBeamData objec function "
                << "is not implemented for parallel run"
                << abort(FatalError);
    }

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    word startPatchName(dict.lookup("startPatchName"));

    polyPatchID startPatch(startPatchName, mesh.boundaryMesh());
    
    if (!startPatch.active())
    {
        FatalErrorIn("writeBeamData::writeBeamData(...)")
          << "Patch name " << startPatchName << " not found."
          << abort(FatalError);
    }
    
    startPatchIndex_ = startPatch.index();
        
    word endPatchName(dict.lookup("endPatchName"));

    polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

    if (!endPatch.active())
    {
        FatalErrorIn("writeBeamData::writeBeamData(...)")
            << "Patch name " << endPatchName << " not found."
            << abort(FatalError);
    }

    endPatchIndex_ = endPatch.index();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::writeBeamData::start()
{
    return false;
}


#if FOAMEXTEND > 40
bool Foam::writeBeamData::execute(const bool forceWrite)
#else
bool Foam::writeBeamData::execute()
#endif
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // const beamModel& beams =
    //     time_.lookupObject<beamModel>("beamProperties");


    if (time_.outputTime())
    {
        OFstream file
        (
            time_.timePath()/"beamShape.dat"
        );

        file.precision(12);

        if (mesh.moving()) // Updated Lagrangian formulatin
        {
            const volVectorField& C = mesh.C();
            const vectorField& CI = C.internalField();

            file << "x" << tab << "y" << tab << "z" << endl; 

            scalar refx = C.boundaryField()[startPatchIndex_][0].x();

            file << C.boundaryField()[startPatchIndex_][0].x() << tab
                 << C.boundaryField()[startPatchIndex_][0].y() << tab
                 << C.boundaryField()[startPatchIndex_][0].z() << endl;

            forAll(CI, cellI)
            {
                file << CI[cellI].x() << tab
                     << CI[cellI].y() << tab
                     << CI[cellI].z() << endl; 
            }

            file << C.boundaryField()[endPatchIndex_][0].x() << tab
                 << C.boundaryField()[endPatchIndex_][0].y() << tab
                 << C.boundaryField()[endPatchIndex_][0].z() << endl;
        }
        else
        {
             const volVectorField& W = 
                 mesh.lookupObject<volVectorField>("W");

            const volVectorField& Theta = 
                mesh.lookupObject<volVectorField>("Theta");
            const vectorField& ThetaI = Theta.internalField();

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

             C += W;
            const vectorField& CI = C.internalField();

           // scalar refx = C.boundaryField()[startPatchIndex_][0].x();
            
            file << "x" << tab << "y" << tab << "z" << endl; 

            file << C.boundaryField()[startPatchIndex_][0].x() << tab
                 << C.boundaryField()[startPatchIndex_][0].y() << tab
                 << C.boundaryField()[startPatchIndex_][0].z() << tab
                 << Theta.boundaryField()[startPatchIndex_][0].x() << tab
                 << Theta.boundaryField()[startPatchIndex_][0].y() << tab
                 << Theta.boundaryField()[startPatchIndex_][0].z() << endl;

            forAll(CI, cellI)
            {
                file << CI[cellI].x() << tab
                     << CI[cellI].y() << tab
                     << CI[cellI].z() << tab
                     << ThetaI[cellI].x() << tab
                     << ThetaI[cellI].y() << tab
                     << ThetaI[cellI].z() << endl; 
            }

            file << C.boundaryField()[endPatchIndex_][0].x() << tab
                 << C.boundaryField()[endPatchIndex_][0].y() << tab
                 << C.boundaryField()[endPatchIndex_][0].z() << tab
                 << Theta.boundaryField()[endPatchIndex_][0].x() << tab
                 << Theta.boundaryField()[endPatchIndex_][0].y() << tab
                 << Theta.boundaryField()[endPatchIndex_][0].z() << endl;
        }

        // Write beam rotational strain
        {
            OFstream file
            (
                time_.timePath()/"beamMaterialRotationalStrains.dat"
            );

            file.precision(12);
            surfaceVectorField Cf = mesh.Cf();

            const surfaceVectorField& K = 
                mesh.lookupObject<surfaceVectorField>("K");

            const surfaceVectorField& KP = 
                mesh.lookupObject<surfaceVectorField>("KP");

            scalar refx = Cf.boundaryField()[startPatchIndex_][0].x();

            file << Cf.boundaryField()[startPatchIndex_][0].x() - refx << tab
                 << -K.boundaryField()[startPatchIndex_][0].x() << tab
                 << -K.boundaryField()[startPatchIndex_][0].y() << tab
                 << -K.boundaryField()[startPatchIndex_][0].z() << tab
                 << -KP.boundaryField()[startPatchIndex_][0].x() << tab
                 << -KP.boundaryField()[startPatchIndex_][0].y() << tab
                 << -KP.boundaryField()[startPatchIndex_][0].z() << endl;

            forAll(K, faceI)
            {
                file << Cf[faceI].x() - refx << tab
                     << K[faceI].x() << tab
                     << K[faceI].y() << tab
                     << K[faceI].z() << tab
                     << KP[faceI].x() << tab
                     << KP[faceI].y() << tab
                     << KP[faceI].z() << endl;
            }
            
            file << Cf.boundaryField()[endPatchIndex_][0].x() - refx << tab
                 << K.boundaryField()[endPatchIndex_][0].x() << tab
                 << K.boundaryField()[endPatchIndex_][0].y() << tab
                 << K.boundaryField()[endPatchIndex_][0].z() << tab
                 << KP.boundaryField()[endPatchIndex_][0].x() << tab
                 << KP.boundaryField()[endPatchIndex_][0].y() << tab
                 << KP.boundaryField()[endPatchIndex_][0].z() << endl;
        }

        // Write beam material moments
        {
            OFstream file
            (
                time_.timePath()/"beamMaterialMoments.txt"
            );

            file.precision(12);
            surfaceVectorField Cf = mesh.Cf();

            const surfaceVectorField& M = 
                mesh.lookupObject<surfaceVectorField>("Mref");

            scalar refx = Cf.boundaryField()[startPatchIndex_][0].x();

            file << Cf.boundaryField()[startPatchIndex_][0].x() - refx << tab
                 << -M.boundaryField()[startPatchIndex_][0].x() << tab
                 << -M.boundaryField()[startPatchIndex_][0].y() << tab
                 << -M.boundaryField()[startPatchIndex_][0].z() << endl;

            forAll(M, faceI)
            {
                file << Cf[faceI].x() - refx << tab
                     << M[faceI].x() << tab
                     << M[faceI].y() << tab
                     << M[faceI].z() << endl;
            }
            
            file << Cf.boundaryField()[endPatchIndex_][0].x() - refx << tab
                 << M.boundaryField()[endPatchIndex_][0].x() << tab
                 << M.boundaryField()[endPatchIndex_][0].y() << tab
                 << M.boundaryField()[endPatchIndex_][0].z() << endl;
        }

        // Write beam axial force
        {
            OFstream file
            (
                time_.timePath()/"beamAxialForce.txt"
            );

            file.precision(12);
            surfaceVectorField Cf = mesh.Cf();

            // const surfaceVectorField& Q = 
            //     mesh.lookupObject<surfaceVectorField>("Q");

            // const surfaceScalarField Qa = mag(Q);

            const surfaceScalarField& Qa = 
                mesh.lookupObject<surfaceScalarField>("Qa");

            scalar refx = Cf.boundaryField()[startPatchIndex_][0].x();

            file << Cf.boundaryField()[startPatchIndex_][0].x() - refx << tab
                 // << Qa.boundaryField()[startPatchIndex_][0] << endl;
                 << -Qa.boundaryField()[startPatchIndex_][0] << endl;

            forAll(Qa, faceI)
            {
                file << Cf[faceI].x() - refx << tab
                     << Qa[faceI] << endl;
            }

            file << Cf.boundaryField()[endPatchIndex_][0].x() - refx << tab
                 << Qa.boundaryField()[endPatchIndex_][0] << endl;
        }
        
        // Write beam curvature
        {
            OFstream file
            (
                time_.timePath()/"beamCurvature.dat"
            );

            file.precision(12);
            const surfaceVectorField& Cf = mesh.Cf();
            const volVectorField& C = mesh.C();

            const volScalarField& curvature = 
                mesh.lookupObject<volScalarField>("curvature");

            scalar refx = Cf.boundaryField()[startPatchIndex_][0].x();

            forAll(curvature, cellI)
            {
                file << C[cellI].x() - refx << tab
                     << curvature[cellI] << endl;
            }
        }
    }

    return false;
}


bool Foam::writeBeamData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
