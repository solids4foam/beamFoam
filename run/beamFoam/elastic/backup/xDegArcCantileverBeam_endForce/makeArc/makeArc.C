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
    
    // Beam mean line displacement
    surfaceVectorField DWf
    (
        IOobject
        (
            "DWf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("DW", dimless, vector::zero)
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

    scalar arcAngle = 8*M_PI;
    
    // Calc displacement and rotation
    boundBox box(mesh.points());
    scalar L = box.max().x() - box.min().x();
    Info << "Beam length: " << L << endl;
    scalar R0 = L/arcAngle;
    Info << "Arc radius: " << R0 << endl;
    
    vectorField& refThetaI = refTheta.internalField();
    vectorField& DWfI = DWf.internalField();
    tensorField& DLambdaI = DLambda.internalField();
    const vectorField& CfI = mesh.Cf().internalField();
    const vectorField& CI = mesh.C().internalField();
    forAll(DWfI, faceI)
    {
        scalar x = CfI[faceI].x();
        scalar alpha = arcAngle*x/L;

        DWfI[faceI] =
            vector(R0*::sin(alpha), R0-R0*::cos(alpha), 0)
          - vector(x, 0, 0);

        DLambdaI[faceI] =
            tensor
            (
                ::cos(alpha), -::sin(alpha), 0,
                ::sin(alpha),  ::cos(alpha), 0,
                0,             0,            1
            );
    }

    forAll(refThetaI, cellI)
    {
        scalar x = CI[cellI].x();
        scalar alpha = arcAngle*x/L;

        refThetaI[cellI] = vector(0, 0, alpha);
    }
    
    forAll(DWf.boundaryField(), patchI)
    {
        vectorField& pRefTheta = refTheta.boundaryField()[patchI];
        vectorField& pDWf = DWf.boundaryField()[patchI];
        tensorField& pDLambda = DLambda.boundaryField()[patchI];
        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];

        forAll(pDWf, faceI)
        {
            scalar x = pCf[faceI].x();
            scalar alpha = arcAngle*x/L;

            pDWf[faceI] =
                vector(R0*::sin(alpha), R0-R0*::cos(alpha), 0)
              - vector(x, 0, 0);

            pDLambda[faceI] =
                tensor
                (
                    ::cos(alpha), -::sin(alpha), 0,
                    ::sin(alpha),  ::cos(alpha), 0,
                    0,             0,            1
                );

            pRefTheta[faceI] = vector(0, 0, alpha);
        }
    }

    refTheta.write();
    DLambda.write();

#   include "updateMeshPoints.H"
  
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
