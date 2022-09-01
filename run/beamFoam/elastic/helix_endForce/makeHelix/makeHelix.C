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
#include "pseudoVector.H"

scalar helixBeta(scalar Lx, scalar R0, scalar L);
scalar f(scalar x, scalar a, scalar b);
tensor R(const vector& n1, const vector& n2);
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

    // scalar arcAngle = M_PI/4;
    // Calc displacement and rotation
    boundBox box(mesh.points());
    scalar L = box.max().x() - box.min().x();
    Info << "Beam length: " << L << endl;

    // Radius of enveloping cylinder
    scalar R0 =
        L
       /(
            6.*::sqrt(sqr(3.*M_PI/4.) + 1.0)
          + (27.*sqr(M_PI)/8.)
           *::log
            (
                (4./(3.*M_PI))
              + ::sqrt(1.0 + sqr(4./(3.*M_PI)))
            )
        );
    Info << "Enveloping cylinder radius: " << R0 << endl;

    vectorField& refTangentI = refTangent.internalField();
    vectorField& DWfI = DWf.internalField();
    tensorField& DLambdaI = DLambda.internalField();
    const vectorField& CfI = mesh.Cf().internalField();
    forAll(DWfI, faceI)
    {
        scalar x = CfI[faceI].x();
        scalar alpha = helixBeta(x, R0, L);

        vector r0
        (
            R0*::sin(alpha),        
            R0*(::cos(alpha) - 1.0),
            R0*(6.0/sqr(9*M_PI))*sqr(alpha)
        );

        DWfI[faceI] = r0 - vector(x, 0, 0);

        tensor T0
        (
            1, 0, 0,
            0,-1, 0,
            0, 0,-1
        );

        vector t
        (
            R0*::cos(alpha),
           -R0*::sin(alpha),
            R0*(12./sqr(9.*M_PI))*alpha
        );
        t /= mag(t);

        vector k(0, 0, 1);
        vector C(0, -R0, 0);
        vector n = C - r0;
        n -= k*(k&n);
        n /= mag(n);
        n -= t*(t&n);
        n /= mag(n);

        // vector n
        // (
        //    -R0*::sin(alpha),
        //    -R0*::cos(alpha),
        //     R0*(12./sqr(9.*M_PI))
        // );
        // n -= t*(t&n);
        // n /= mag(n);

        vector bn = (t^n);

        tensor T1
        (
            t.x(), n.x(), bn.x(),
            t.y(), n.y(), bn.y(),
            t.z(), n.z(), bn.z()
        );

        DLambdaI[faceI] = R(T0, T1); //R(T, t); //(t*T);
        refTangentI[faceI] = t;

        // vector newT = (DLambdaI[faceI] & vector(1, 0, 0));
        // vector newCor = (DLambdaI[faceI] & vector(1, 0, 0));
        // vector newT0 = (T1 & newCor);
        // vector newT = (T1 & newT0);
        // vector newT = newCor.x()*t + newCor.y()*n + newCor.z()*bn;
        // Info << t << ", " << newCor << ", " << newT << endl;
        // Info << "AAT: " << (DLambdaI[faceI] & DLambdaI[faceI].T()) << endl;
    }

    forAll(DWf.boundaryField(), patchI)
    {
        vectorField& pRefTangent = refTangent.boundaryField()[patchI];
        vectorField& pDWf = DWf.boundaryField()[patchI];
        tensorField& pDLambda = DLambda.boundaryField()[patchI];
        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];

        forAll(pDWf, faceI)
        {
            scalar x = pCf[faceI].x();
            scalar alpha = helixBeta(x, R0, L);

            vector r0
            (
                R0*::sin(alpha),        
                R0*(::cos(alpha) - 1.0),
                R0*(6.0/sqr(9*M_PI))*sqr(alpha)
            );

            pDWf[faceI] = r0 - vector(x, 0, 0);

            tensor T0
            (
                1, 0, 0,
                0,-1, 0,
                0, 0,-1
            );

            vector t
            (
                R0*::cos(alpha),
               -R0*::sin(alpha),
                R0*(12./sqr(9.*M_PI))*alpha
            );
            t /= mag(t);

            vector k(0, 0, 1);
            vector C(0, -R0, 0);
            vector n = C - r0;
            n -= k*(k&n);
            n /= mag(n);
            n -= t*(t&n);
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
    }

    DLambda.write();
    refTangent.write();

#   include "updateMeshPoints.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    // {
    //   vector n1(1, 0, 0);
    //   vector n2(-1, 0, 0);
    //   tensor rot = R(n1, n2);

    //   Info << rot << endl;
    // }
    
    return(0);
}

scalar helixBeta(scalar Lx, scalar R0, scalar L)
{
    scalar a = sqr(R0);
    scalar b = sqr(12*R0/(81*sqr(M_PI)));

    scalar betaMax = 9*M_PI;

    // scalar len = f(betaMax, a, b) - f(0, a, b);
    // Info << "len = " << len << endl;

    if (Lx < SMALL)
    {
        return 0;
    }

    if (Lx > (L-SMALL))
    {
        return betaMax;
    }

    scalar residual = GREAT;
    scalar x0 = 0;
    scalar x1 = betaMax;
    scalar x2 = 0;
    do
    {
        scalar f0 = f(x0, a, b) - f(0, a, b) - Lx;
        scalar f1 = f(x1, a, b) - f(0, a, b) - Lx;

        x2 = x1 - f1*(x1-x0)/(f1-f0);
        
        residual = mag(x2-x1)/betaMax;

        x0 = x1;
        x1 = x2;
    }
    while(residual > 1e-4);

    // scalar len = f(x2, a, b) - f(0, a, b);
    // Info << Lx << ", " << len << endl;
    
    return x2;
}

scalar f(scalar x, scalar a, scalar b)
{
    scalar result =
      0.5
     *(
          x*::sqrt(a+b*sqr(x))
        + a*::log(::sqrt(b)*::sqrt(a+b*sqr(x)) + b*x)
         /::sqrt(b)
      );
 
    return result;
}

inline tensor R
(
    const vector& n1,
    const vector& n2
)
{
  // Info << magSqr(n1 ^ n2) << endl;
  
    tensor result =
        (n1 & n2)*I
      + (1 - (n1 & n2))*sqr(n1 ^ n2)/(magSqr(n1 ^ n2) + VSMALL)
      + (n2*n1 - n1*n2);

    // Info << det(result) << endl;
    
    return result;
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
