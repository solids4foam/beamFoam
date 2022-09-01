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

Application
    setSecondBeam

Description
    Calculate initial position and shape of second beam

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"
#include "boundBox.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    surfaceVectorField  refTangent
    (
        IOobject
        (
            "refTangent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("x", dimless, vector(1, 0, 0))
    );

    // Beam mean line displacement
    surfaceVectorField refWf
    (
        IOobject
        (
            "refWf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // Beam mean line displacement
    volVectorField refW
    (
        IOobject
        (
            "refW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );
    
    // Beam cross-section rotation
    surfaceTensorField refLambda
    (
        IOobject
        (
            "refLambda",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );

    volTensorField refRM
    (
        IOobject
        (
            "refRM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );
    tensorField& refRMI = refRM.internalField();

    // Set rotation and translation
    scalar alphaY = 90*M_PI/180;

    tensor Ry =
        tensor
        (
            ::cos(alphaY), 0, ::sin(alphaY),
            0, 1, 0,
           -::sin(alphaY), 0, ::cos(alphaY)
        );

    // vector DR(0, 0, 0);
    vector DR(0.25, 0.01, 0.25);
    
    // refLambda = Ry;
    // refRM = Ry;

    vectorField& refWfI = refWf.internalField();
    vectorField& refWI = refW.internalField();
    vectorField& refTangentI = refTangent.internalField();
    tensorField& refLambdaI = refLambda.internalField();

    const vectorField& CfI = mesh.Cf().internalField();
    const vectorField& CI = mesh.C().internalField();

    vectorField newCfI = (Ry & CfI);

    const labelList& nei = mesh.neighbour();
    
    forAll(refWfI, faceI)
    {
        label I = mesh.cellZones().whichZone(nei[faceI]);

        // Only second beam
        if (I == 1)
        {
            refWfI[faceI] = (newCfI[faceI] + DR - CfI[faceI]);

            refLambdaI[faceI] = Ry;
            
            refTangentI[faceI] = (refLambdaI[faceI] & refTangentI[faceI]);
        }
    }


    vectorField newCI = (Ry & CI);
    
    forAll(refRMI, cellI)
    {
        label I = mesh.cellZones().whichZone(cellI);

        // Only second beam
        if (I == 1)
        {
            refWI[cellI] = (newCI[cellI] + DR - CI[cellI]);
            refRMI[cellI] = Ry;
        }
    }
    
    forAll(refWf.boundaryField(), patchI)
    {
        vectorField& pRefWf = refWf.boundaryField()[patchI];
        vectorField& pRefW = refW.boundaryField()[patchI];
        vectorField& pRefTangent = refTangent.boundaryField()[patchI];
        tensorField& pRefLambda = refLambda.boundaryField()[patchI];
        tensorField& pRefRM = refRM.boundaryField()[patchI];

        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];

        const vectorField newpCf = (Ry & pCf);

        const labelList faceCells =
            mesh.boundary()[patchI].faceCells();
        
        forAll(pRefWf, faceI)
        {
            label I = mesh.cellZones().whichZone(faceCells[faceI]);

            if (I == 1)
            {
                pRefWf[faceI] = (newpCf[faceI] + DR - pCf[faceI]);
                pRefW[faceI] = pRefWf[faceI];

                pRefLambda[faceI] = Ry;
                pRefRM[faceI] = Ry;

                pRefTangent[faceI] = (pRefLambda[faceI] & pRefTangent[faceI]);
            }
        }
    }

    refWf.write();
    refW.write();
    refLambda.write();
    refTangent.write();
    refRM.write();

#   include "updateMeshPoints.H"
  
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
