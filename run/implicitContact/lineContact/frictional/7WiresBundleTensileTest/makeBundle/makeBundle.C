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
    makeBundle

Description
    Deforming streight beams to bundle

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"
#include "boundBox.H"
#include "pseudoVector.H"

tensor R
(
    const tensor& T0,
    const tensor& T1
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    // Beam cross-section rotation
    surfaceVectorField  refTangent
    (
        IOobject
        (
            "refTangent",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector(1, 0, 0))
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

    // Calc displacement and rotation
    boundBox box(mesh.points());
    scalar L = box.max().x() - box.min().x();
    Info << "Beam length: " << L << endl;

    // Radius of enveloping cylinder
    scalar R0 = 0.02;
    Info << "Enveloping cylinder radius: " << R0 << endl;

    // Slope
    scalar h = ::sqrt(sqr(L/(8*M_PI)) - sqr(R0));
    Info << "Helix slope: " << h << endl;

    vectorField& refTangentI = refTangent.internalField();
    vectorField& DWfI = DWf.internalField();
    tensorField& DLambdaI = DLambda.internalField();
    const vectorField& CfI = mesh.Cf().internalField();

    const labelList& nei = mesh.neighbour();
    // const cellZone& cz = mesh.cellZones()[bI];

    forAll(DWfI, faceI)
    {
        label I = mesh.cellZones().whichZone(nei[faceI]);
        scalar delta = 2*M_PI/6;

        if (I>0)
        {
            scalar x = CfI[faceI].x();
            scalar phi = x/::sqrt(sqr(h) + sqr(R0));

            vector r0
            (
                h*phi,
                R0*::cos(phi + (I-1)*delta),
                R0*::sin(phi + (I-1)*delta)
            );

            DWfI[faceI] = r0 - CfI[faceI];

            tensor RT
            (
                1, 0,                    0,
                0, ::cos(-(I-1)*delta), -::sin(-(I-1)*delta),
                0, ::sin(-(I-1)*delta),  ::cos(-(I-1)*delta)
            );

            vector i(1,  0, 0);
            vector j(0, -1, 0);
            j = (RT & j);
            vector k(0, 0, -1);
            k = (RT & k);
            
            tensor T0
            (
                i.x(), j.x(), k.x(),
                i.y(), j.y(), k.y(),
                i.z(), j.z(), k.z()
            );

            vector t
            (
                h,
               -R0*::sin(phi + (I-1)*delta),
                R0*::cos(phi + (I-1)*delta)
            );
            t /= mag(t);

            vector n
            (
                0,
               -R0*::cos(phi + (I-1)*delta),
               -R0*::sin(phi + (I-1)*delta)
            );
            n /= mag(n);

            vector bn = (t^n);

            tensor T1
            (
                t.x(), n.x(), bn.x(),
                t.y(), n.y(), bn.y(),
                t.z(), n.z(), bn.z()
            );

            DLambdaI[faceI] = R(T0, T1);
            refTangentI[faceI] = t;
        }
        else
        {
            scalar x = CfI[faceI].x();
            scalar phi = x/::sqrt(sqr(h) + sqr(R0));

            vector r0
            (
                h*phi,
                R0*::cos(phi),
                R0*::sin(phi)
            );

            DWfI[faceI] = vector(r0.x(), 0, 0) - CfI[faceI];
        }
    }

    forAll(DWf.boundaryField(), patchI)
    {
        vectorField& pRefTangent = refTangent.boundaryField()[patchI];
        vectorField& pDWf = DWf.boundaryField()[patchI];
        tensorField& pDLambda = DLambda.boundaryField()[patchI];
        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];

        const labelList faceCells =
            mesh.boundary()[patchI].faceCells();
        
        forAll(pDWf, faceI)
        {
            label I = mesh.cellZones().whichZone(faceCells[faceI]);
            scalar delta = 2*M_PI/6;
        
            if (I>0)
            {
                scalar x = pCf[faceI].x();
                scalar phi = x/::sqrt(sqr(h) + sqr(R0));

                vector r0
                (
                    h*phi,
                    R0*::cos(phi + (I-1)*delta),
                    R0*::sin(phi + (I-1)*delta)
                );

                pDWf[faceI] = r0 - pCf[faceI];

                tensor RT
                (
                    1, 0,                    0,
                    0, ::cos(-(I-1)*delta), -::sin(-(I-1)*delta),
                    0, ::sin(-(I-1)*delta),  ::cos(-(I-1)*delta)
                );

                vector i(1,  0, 0);
                vector j(0, -1, 0);
                j = (RT & j);
                vector k(0, 0, -1);
                k = (RT & k);

                tensor T0
                (
                    i.x(), j.x(), k.x(),
                    i.y(), j.y(), k.y(),
                    i.z(), j.z(), k.z()
                );
            
                vector t
                (
                    h,
                   -R0*::sin(phi + (I-1)*delta),
                    R0*::cos(phi + (I-1)*delta)
                );
                t /= mag(t);

                vector n
                (
                    0,
                   -R0*::cos(phi + (I-1)*delta),
                   -R0*::sin(phi + (I-1)*delta)
                );
                n /= mag(n);

                vector bn = (t^n);

                tensor T1
                (
                    t.x(), n.x(), bn.x(),
                    t.y(), n.y(), bn.y(),
                    t.z(), n.z(), bn.z()
                );
        
                pDLambda[faceI] = R(T0, T1);
                pRefTangent[faceI] = t;
            }
            else
            {
                scalar x = pCf[faceI].x();
                scalar phi = x/::sqrt(sqr(h) + sqr(R0));

                vector r0
                (
                    h*phi,
                    R0*::cos(phi),
                    R0*::sin(phi)
                );

                pDWf[faceI] = vector(r0.x(), 0, 0) - pCf[faceI];
            }
        }
    }

    DLambda.write();
    refTangent.write();

#   include "updateMeshPoints.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}

inline tensor R
(
    const tensor& T0,
    const tensor& T1
)
{
    vector X0(T0.xx(), T0.yx(), T0.zx());
    vector Y0(T0.xy(), T0.yy(), T0.zy());
    vector Z0(T0.xz(), T0.yz(), T0.zz());

    vector X1(T1.xx(), T1.yx(), T1.zx());
    vector Y1(T1.xy(), T1.yy(), T1.zy());
    vector Z1(T1.xz(), T1.yz(), T1.zz());

    // tensor result
    // (
    //    (X0&X1), (X0&Y1), (X0&Z0), 
    //    (Y0&X1), (Y0&Y1), (Y0&Z0), 
    //    (Z0&X1), (Z0&Y1), (Z0&Z0)      
    // );
    tensor result
    (
       (X1&X0), (X1&Y0), (X1&Z0), 
       (Y1&X0), (Y1&Y0), (Y1&Z0), 
       (Z1&X0), (Z1&Y0), (Z1&Z0)      
    );

    result = (T1 & (T1 & result));
    
    return result;
}

// ************************************************************************* //
