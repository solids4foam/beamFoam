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
    makeCircularArc

Description
    Deforming straight beam geometry to the 215deg arc shape

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

    // Beam cross-section rotation
    volVectorField  refTheta
    (
        IOobject
        (
            "refTheta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    );

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

    // Beam cross-section rotation
    surfaceTensorField DLambda
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

    //argList::validOptions.insert("arcAngleInDegrees", "scalar");
    
    //scalar arcAngle = (readScalar(args.optionLookup("arcAngleInDegrees")()))*M_PI/180;
    
    scalar arcAngle = 215*M_PI/180;
    
    // Calc displacement and rotation
    boundBox box(mesh.points());
    scalar L = box.max().x() - box.min().x();
    Info << "Beam length: " << L << endl;
    scalar R0 = L/arcAngle;
    Info << "Arc radius: " << R0 << endl;
    
    vectorField& refThetaI = refTheta.internalField();
    vectorField& refWfI = refWf.internalField();
    tensorField& DLambdaI = DLambda.internalField();
    vectorField& refTangentI = refTangent.internalField();
    const vectorField& CfI = mesh.Cf().internalField();
    const vectorField& CI = mesh.C().internalField();
    
    forAll(refWfI, faceI)
    {
        scalar x = CfI[faceI].x();
        scalar beta = arcAngle*x/L;
        scalar alpha = 0.5*(2*M_PI - arcAngle);
        scalar phi = M_PI - (alpha + beta);

        refWfI[faceI] =
            vector
            (
            	R0*(1 - (::cos(beta) + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::sin(alpha),
            	R0*(1 - (::cos(beta) + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::cos(alpha) + R0*::sin(beta)/::sin(alpha),
            	0
            )
          - vector(x, 0, 0);

        DLambdaI[faceI] =
            tensor
            (
                ::cos(phi), -::sin(phi), 0,
                ::sin(phi),  ::cos(phi), 0,
                0,             0,            1
            );

        refTangentI[faceI] = (DLambdaI[faceI] & refTangentI[faceI]);
    }
    
    
    forAll(CI, cellI)
    {
        scalar x = CI[cellI].x();
        scalar beta = arcAngle*x/L;
        scalar alpha = 0.5*(2*M_PI - arcAngle);
        scalar phi = M_PI - (alpha + beta);

        refRMI[cellI] =
            tensor
            (
                ::cos(phi), -::sin(phi), 0,
                ::sin(phi),  ::cos(phi), 0,
                0,             0,            1
            );
    }


    forAll(refThetaI, cellI)
    {
        scalar x = CI[cellI].x();
        scalar beta = arcAngle*x/L;

        refThetaI[cellI] = vector(0, 0, beta);
    }
    
    forAll(refWf.boundaryField(), patchI)
    {
        vectorField& pRefTheta = refTheta.boundaryField()[patchI];
        vectorField& pRefWf = refWf.boundaryField()[patchI];
        tensorField& pDLambda = DLambda.boundaryField()[patchI];
        vectorField& pRefTangent = refTangent.boundaryField()[patchI];
        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];
        
        tensorField& pRefRM = refRM.boundaryField()[patchI];

        forAll(pRefWf, faceI)
        {
            scalar x = pCf[faceI].x();
            scalar beta = arcAngle*x/L;
            scalar alpha = 0.5*(2*M_PI - arcAngle);
            scalar phi = M_PI - (alpha + beta);

            pRefWf[faceI] =
             	vector
            (
            	R0*(1 - (::cos(beta) + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::sin(alpha),
            	R0*(1 - (::cos(beta) + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::cos(alpha) + R0*::sin(beta)/::sin(alpha),
            	0
            )
              - vector(x, 0, 0);

            pDLambda[faceI] =
                tensor
                (
                    ::cos(phi), -::sin(phi), 0,
                    ::sin(phi),  ::cos(phi), 0,
                    0,             0,        1
                );

            pRefTheta[faceI] = vector(0, 0, beta);

            pRefTangent[faceI] = (pDLambda[faceI] & pRefTangent[faceI]);
            
            pRefRM[faceI] = pDLambda[faceI];
        }
    }
	
	refWf.write();
    DLambda.write();
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
