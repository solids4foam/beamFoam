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
    makeArc

Description
    Deforming streight beam geometry to the 45deg arc shape

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

    // // Beam cross-section rotation
    // volVectorField  refTheta
    // (
    //     IOobject
    //     (
    //         "refTheta",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedVector("zero", dimless, vector::zero)
    // );

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

    scalar arcAngle = -M_PI/2;

    // Calc displacement and rotation
    boundBox box(mesh.points());
    scalar L = box.max().x() - box.min().x();
    Info << "Beam length: " << L << endl;
    scalar R0 = L/arcAngle;
    Info << "Arc radius: " << R0 << endl;

    // vectorField& refThetaI = refTheta.internalField();
    vectorField& refWfI = refWf.internalField();
    vectorField& refWI = refW.internalField();
    tensorField& refLambdaI = refLambda.internalField();
    vectorField& refTangentI = refTangent.internalField();
    const vectorField& CfI = mesh.Cf().internalField();
    const vectorField& CI = mesh.C().internalField();

    forAll(refWfI, faceI)
    {
        scalar x = CfI[faceI].x();
        scalar alpha = arcAngle*x/L;

        refWfI[faceI] =
            vector(R0*::sin(alpha), R0-R0*::cos(alpha), 0)
          - vector(x, 0, 0);

        refLambdaI[faceI] =
            tensor
            (
                ::cos(alpha), -::sin(alpha), 0,
                ::sin(alpha),  ::cos(alpha), 0,
                0,             0,            1
            );

        refTangentI[faceI] = (refLambdaI[faceI] & refTangentI[faceI]);
    }

    forAll(CI, cellI)
    {
        scalar x = CI[cellI].x();
        scalar alpha = arcAngle*x/L;

        refWI[cellI] =
            vector(R0*::sin(alpha), R0-R0*::cos(alpha), 0)
          - vector(x, 0, 0);
        
        refRMI[cellI] =
            tensor
            (
                ::cos(alpha), -::sin(alpha), 0,
                ::sin(alpha),  ::cos(alpha), 0,
                0,             0,            1
            );
    }

    // forAll(refThetaI, cellI)
    // {
    //     scalar x = CI[cellI].x();
    //     scalar alpha = arcAngle*x/L;

    //     refThetaI[cellI] = vector(0, 0, alpha);
    // }
    
    forAll(refWf.boundaryField(), patchI)
    {
        // vectorField& pRefTheta = refTheta.boundaryField()[patchI];
        vectorField& pRefWf = refWf.boundaryField()[patchI];
        vectorField& pRefW = refW.boundaryField()[patchI];
        tensorField& pRefLambda = refLambda.boundaryField()[patchI];
        vectorField& pRefTangent = refTangent.boundaryField()[patchI];
        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];

        tensorField& pRefRM = refRM.boundaryField()[patchI];

        forAll(pRefWf, faceI)
        {
            scalar x = pCf[faceI].x();
            scalar alpha = arcAngle*x/L;

            pRefWf[faceI] =
                vector(R0*::sin(alpha), R0-R0*::cos(alpha), 0)
              - vector(x, 0, 0);

            pRefW[faceI] = pRefWf[faceI];
                
            pRefLambda[faceI] =
                tensor
                (
                    ::cos(alpha), -::sin(alpha), 0,
                    ::sin(alpha),  ::cos(alpha), 0,
                    0,             0,            1
                );

            // pRefTheta[faceI] = vector(0, 0, alpha);

            pRefTangent[faceI] = (pRefLambda[faceI] & pRefTangent[faceI]);

            //
            pRefRM[faceI] = pRefLambda[faceI];
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
