/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description

Author
    Zeljko Tukovic, FSB Zagreb

\*---------------------------------------------------------------------------*/

#include "pseudoVector.H"
#include "emptyFvPatchFields.H"
#include "spinTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<volVectorField> pseudoVector(const volTensorField& rotationMatrix)
{
    tmp<volVectorField> tresult
    (    
        new volVectorField
        (
            IOobject
            (
                "pseudoVector("+rotationMatrix.name()+")",
                rotationMatrix.time().timeName(),
                rotationMatrix.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotationMatrix.mesh(),
            dimensionedVector("0", dimless, vector::zero)
        )
    );

    vectorField& resultI = tresult().internalField();

    const tensorField& rotationMatrixI = rotationMatrix.internalField();
    
    forAll(resultI, cellI)
    {
        resultI[cellI] = pseudoVector(rotationMatrixI[cellI]);
    }

    forAll (rotationMatrix.boundaryField(), patchI)
    {
        vectorField& pResultI = tresult().boundaryField()[patchI];
        
        const tensorField& pRotationMatrixI =
            rotationMatrix.boundaryField()[patchI];
        
        forAll(pResultI, faceI)
        {
            pResultI[faceI] = pseudoVector(pRotationMatrixI[faceI]);
        }
    }

    return tresult;
}
    
tmp<vectorField> pseudoVector(const tensorField& rotationMatrix)
{
    tmp<vectorField> tresult
    (
        new vectorField(rotationMatrix.size(), vector::zero)
    );

    vectorField& result = tresult();

    forAll(result, cellI)
    {
        result[cellI] = pseudoVector(rotationMatrix[cellI]);
    }

    return tresult;
}
  
vector pseudoVector(const tensor& R)
{
    vector result = vector::zero;

    if (false)
    {
        vector rotAxis
        (
            R.zy() - R.yz(),
            R.xz() - R.zx(),
            R.yx() - R.xy()
        );

        scalar magRotAngle = asin(mag(rotAxis)/2);

        rotAxis /= mag(rotAxis) + SMALL;

        result = magRotAngle*rotAxis;

        return result;
    }

    scalar q0 = 0;
    vector q = vector::zero;

    scalarField D(4, 0);
    D[0] = R.xx();
    D[1] = R.yy();
    D[2] = R.zz();
    D[3] = tr(R);

    scalar a = max(D);

    // if ( mag(a-D[3]) < SMALL )
    // {
    //     q0 = 0.5*::sqrt(1+a);
    //     q.x() = (R.zy()-R.yz())/(4*q0);
    //     q.y() = (R.xz()-R.zx())/(4*q0);
    //     q.z() = (R.yz()-R.zy())/(4*q0);
    // }
    // else
    {
        label i = -1;
        for (label j=0; j<4; j++)
        {
            if ( mag(a-D[j]) < SMALL )
            {
                i = j;
                break;
            }
        }

        // Info << D << endl;
        // Info << "a = " << a << endl;
 
        switch(i)
        {
            case 0:
                q.x() = sqrt(a/2 + (1-tr(R))/4);
                q0 = (R.zy()-R.yz())/(4*q.x());
                q.y() = (R.yx()+R.xy())/(4*q.x());
                q.z() = (R.zx()+R.xz())/(4*q.x());
                break;
            case 1:
                q.y() = sqrt(a/2 + (1-tr(R))/4);
                q0 = (R.xz()-R.zx())/(4*q.y());
                q.z() = (R.zy()+R.yz())/(4*q.y());
                q.x() = (R.xy()+R.yx())/(4*q.y());
                break;
            case 2:
                q.z() = sqrt(a/2 + (1-tr(R))/4);
                q0 = (R.yz()-R.zy())/(4*q.z());
                q.x() = (R.xz()+R.zx())/(4*q.z());
                q.y() = (R.yz()+R.zy())/(4*q.z());
                break;
            case 3:
                q0 = 0.5*::sqrt(1+a);
                q.x() = (R.zy()-R.yz())/(4*q0);
                q.y() = (R.xz()-R.zx())/(4*q0);
                q.z() = (R.yx()-R.xy())/(4*q0);
                break;
            default:
                FatalErrorIn
                (
                    "vector pseudoVector(const tensor& R)"
                )
                  << "Error in calculation of quaternions"
                  << "for rotation matrix " << R
                  << abort(FatalError);
        }
    }

    scalar theta = atan(mag(q)/q0)*2;

    if (mag(theta) > SMALL)
    {
        result = theta*(q/mag(q));
    }

    return result;
}

    
tmp<tensorField> rotationMatrix(const vectorField& rotationAngle)
{
    tmp<tensorField> tRotMat
    (    
        new tensorField(rotationAngle.size(), tensor::I)
    );

    if (false)
    {
        scalarField magAngle = mag(rotationAngle);
        vectorField rotAxis = rotationAngle/(magAngle + SMALL);
        tensorField rotAxisHat = spinTensor(rotAxis);

        tRotMat() =
            Foam::cos(magAngle)*I +
            Foam::sin(magAngle)*rotAxisHat +
            (1.0 - Foam::cos(magAngle))*(rotAxis*rotAxis);
    }
    else
    {
        scalarField magAngle = mag(rotationAngle) + SMALL;
        tensorField angleHat = spinTensor(rotationAngle);
    
        tRotMat() =
            I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
           *(angleHat & angleHat);
    }
    
    return tRotMat;
}

    
tmp<volTensorField> rotationMatrix(const volVectorField& rotationAngle)
{
    tmp<volTensorField> tRotMat
    (    
        new volTensorField
        (
            IOobject
            (
                "RotMat("+rotationAngle.name()+")",
                rotationAngle.time().timeName(),
                rotationAngle.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotationAngle.mesh(),
            dimensionedTensor("0", dimless, tensor::I)
        )
    );

    if (false)
    {
        volScalarField magAngle = mag(rotationAngle);
        volVectorField rotAxis = rotationAngle/(magAngle + SMALL);
        volTensorField rotAxisHat = spinTensor(rotAxis);

        tRotMat() =
            Foam::cos(magAngle)*I +
            Foam::sin(magAngle)*rotAxisHat +
            (1.0 - Foam::cos(magAngle))*(rotAxis*rotAxis);
    }
    else
    {
        volScalarField magAngle = mag(rotationAngle) + SMALL;
        volTensorField angleHat = spinTensor(rotationAngle);
    
        tRotMat() =
            I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
           *(angleHat & angleHat);
    }
    
    return tRotMat;
}

  
tmp<surfaceTensorField> rotationMatrix(const surfaceVectorField& rotationAngle)
{
    tmp<surfaceTensorField> tRotMat
    (    
        new surfaceTensorField
        (
            IOobject
            (
                "RotMat("+rotationAngle.name()+")",
                rotationAngle.time().timeName(),
                rotationAngle.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotationAngle.mesh(),
            dimensionedTensor("0", dimless, tensor::I)
        )
    );

    if (false)
    {
        surfaceScalarField magAngle = mag(rotationAngle);
        surfaceVectorField rotAxis = rotationAngle/(magAngle + SMALL);
        surfaceTensorField rotAxisHat = spinTensor(rotAxis);

        tRotMat() =
            Foam::cos(magAngle)*I +
            Foam::sin(magAngle)*rotAxisHat +
            (1.0 - Foam::cos(magAngle))*(rotAxis*rotAxis);
    }
    else
    {
        surfaceScalarField magAngle = mag(rotationAngle) + SMALL;
        surfaceTensorField angleHat = spinTensor(rotationAngle);
    
        tRotMat() =
            I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
           *(angleHat & angleHat);
    }
    return tRotMat;
}
    
tensor rotationMatrix(const vector& rotationAngle)
{
    tensor R = tensor::zero;
    
    if (false)
    {
        scalar magAngle = mag(rotationAngle);
        vector rotAxis = rotationAngle/(magAngle + SMALL);
        tensor rotAxisHat = spinTensor(rotAxis);

        // Rodrigues formula
        R =
            Foam::cos(magAngle)*I +
            Foam::sin(magAngle)*rotAxisHat +
            (1.0 - Foam::cos(magAngle))*(rotAxis*rotAxis);
    }
    else
    {
        scalar magAngle = mag(rotationAngle) + SMALL;
        tensor angleHat = spinTensor(rotationAngle);

        // Rodrigues formula
        R =
            I + (Foam::sin(magAngle)/magAngle)*angleHat
          + ((1.0-Foam::cos(magAngle))/sqr(magAngle))
           *(angleHat & angleHat);
    }
    
    return R;
}

tmp<surfaceTensorField> interpolateRotationMatrix
(
    const volTensorField& R
)
{
    tmp<surfaceTensorField> tRf
    (    
        new surfaceTensorField
        (
            IOobject
            (
                "Rf",
                R.time().timeName(),
                R.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R.mesh(),
            dimensionedTensor("0", dimless, tensor::I)
        )
    );

    surfaceTensorField& Rf = tRf();

    forAll(Rf.boundaryField(), patchI)
    {
    	// For boundary faces, why are we substituting the cell centre rotation values?
        Rf.boundaryField()[patchI] =
            R.boundaryField()[patchI];
    }

    const tensorField& RI = R.internalField();

    tensorField& RfI = Rf.internalField();
    
    const fvMesh& mesh = R.mesh();
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();

    forAll(RfI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);
        
        vector DRotAngle = pseudoVector(DR);

        // scalar s = 0.5/deltaCoeffsI[faceI];
        RfI[faceI] = RI[own[faceI]] & rotationMatrix(0.5*DRotAngle);
      
    }

    return tRf;
}


tmp<surfaceVectorField> meanLineCurvature
(
    const volTensorField& R
)
{
    tmp<surfaceVectorField> tKf
    (    
        new surfaceVectorField
        (
            IOobject
            (
                "Kf",
                R.time().timeName(),
                R.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            R.mesh(),
            dimensionedVector("0", dimless/dimLength, vector::zero)
        )
    );

    surfaceVectorField& Kf = tKf();

    const tensorField& RI = R.internalField();

    vectorField& KfI = Kf.internalField();
    
    const fvMesh& mesh = R.mesh();
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();
    const scalarField& deltaCoeffsI = mesh.deltaCoeffs();

    forAll(KfI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);

        vector DRotAngle = pseudoVector(DR);

        KfI[faceI] = DRotAngle*deltaCoeffsI[faceI];
    }

    forAll(Kf.boundaryField(), patchI)
    {
        vectorField& pKf = Kf.boundaryField()[patchI];
        const tensorField& pR = R.boundaryField()[patchI];

        labelList faceCells = mesh.boundary()[patchI].faceCells();
        const scalarField& pDeltaCoeffs =
            mesh.deltaCoeffs().boundaryField()[patchI];

        forAll(pKf, faceI)
        {
            tensor DR = (RI[faceCells[faceI]].T() & pR[faceI]);
            vector DRotAngle = pseudoVector(DR);

            tensor RN = (RI[faceCells[faceI]] & rotationMatrix(2*DRotAngle));
            
            DR = (RI[faceCells[faceI]].T() & RN);
            DRotAngle = pseudoVector(DR);
                        
            pKf[faceI] = DRotAngle*pDeltaCoeffs[faceI]/2;
        }
    }

    return tKf;
}

// This member function is not being used anymore to interpolate the rotation matrix
void interpolateRotationMatrix
(
    const volTensorField& R,
    surfaceTensorField& Rf,
    surfaceVectorField& snGradRotAngle
)
{
    forAll(Rf.boundaryField(), patchI)
    {
        Rf.boundaryField()[patchI] = R.boundaryField()[patchI];
    }

    const tensorField& RI = R.internalField();

    tensorField& RfI = Rf.internalField();
    vectorField& snGradRotAngleI = snGradRotAngle;

    const fvMesh& mesh = R.mesh();
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();
    const scalarField& deltaCoeffsI = mesh.deltaCoeffs();

    forAll(snGradRotAngleI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);

        vector DRotAngle = pseudoVector(DR);

        snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];

        // scalar s = 0.5/deltaCoeffsI[faceI];
        RfI[faceI] = RI[own[faceI]] & rotationMatrix(0.5*DRotAngle);
    }

    forAll(snGradRotAngle.boundaryField(), patchI)
    {
        const tensorField& RI = R.boundaryField()[patchI];

        // tensorField& RfI = Rf.boundaryField()[patchI];
        vectorField& snGradRotAngleI =
            snGradRotAngle.boundaryField()[patchI];

        const labelList& own =
            mesh.boundary()[patchI].faceCells();
        const scalarField& deltaCoeffsI =
            mesh.deltaCoeffs().boundaryField()[patchI];

        forAll(snGradRotAngleI, faceI)
        {
            tensor DR = (R[own[faceI]].T() & RI[faceI]);

            vector DRotAngle = pseudoVector(DR);

            snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];
        }
    }
}

// This member function is not being used anymore to interpolate the rotation matrix
void interpolateRotationMatrix
(
    const beamModel& bm,
    const surfaceTensorField& Rf,
    volTensorField& R
)
{
    forAll(R.boundaryField(), patchI)
    {
        if (!R.boundaryField()[patchI].coupled())
        {
            R.boundaryField()[patchI] = Rf.boundaryField()[patchI];
        }
    }

    const tensorField& RfI = Rf.internalField();

    const labelList& upperNeiFaces = bm.upperNeiCellFaces();
    const labelList& lowerNeiFaces = bm.lowerNeiCellFaces();

    // Pout << "y1" << endl;
    // sleep(5);

    tensorField& RI = R.internalField();

    const fvMesh& mesh = R.mesh();
    // const unallocLabelList& own = mesh.owner();
    // const unallocLabelList& nei = mesh.neighbour();

    forAll(RI, cellI)
    {
        if (lowerNeiFaces[cellI] < 0)
        {
            label startPatchIndex = -lowerNeiFaces[cellI] - 1;

            // if (Rf.boundaryField()[startPatchIndex].size() > 1)
            // {
            //     FatalErrorIn
            //     (
            //         "interpolateRotationMatrix(...)"
            //     )
            //       << "Face to cell rotation matrix interpolation "
            //       << "is not implemented for multi-beam problem"
            //       << abort(FatalError);
            // }

            label bFaceI =
                findIndex
                (
                    mesh.boundary()[startPatchIndex].faceCells(),
                    cellI
                );

            tensor DR =
            (
                Rf.boundaryField()[startPatchIndex][bFaceI].T()
              & RfI[upperNeiFaces[cellI]]
            );

            vector DRotAngle = pseudoVector(DR);

            RI[cellI] =
            (
                Rf.boundaryField()[startPatchIndex][0].T()
              & rotationMatrix(0.5*DRotAngle)
            );
        }
        else if (upperNeiFaces[cellI] < 0)
        {
            label endPatchIndex = -upperNeiFaces[cellI] - 1;

            // if (Rf.boundaryField()[endPatchIndex].size() > 1)
            // {
            //     FatalErrorIn
            //     (
            //         "interpolateRotationMatrix(...)"
            //     )
            //       << "Face to cell rotation matrix interpolation "
            //       << "is not implemented for mutiple-beam problem"
            //       << abort(FatalError);
            // }

            label bFaceI =
                findIndex
                (
                    mesh.boundary()[endPatchIndex].faceCells(),
                    cellI
                );

            tensor DR =
            (
                RfI[lowerNeiFaces[cellI]].T()
              & Rf.boundaryField()[endPatchIndex][bFaceI]
            );

            vector DRotAngle = pseudoVector(DR);

            RI[cellI] =
            (
                RfI[lowerNeiFaces[cellI]]
              & rotationMatrix(0.5*DRotAngle)
            );
        }
        else
        {
            tensor DR =
            (
                RfI[lowerNeiFaces[cellI]].T()
              & RfI[upperNeiFaces[cellI]]
            );

            vector DRotAngle = pseudoVector(DR);

            RI[cellI] =
            (
                RfI[lowerNeiFaces[cellI]]
              & rotationMatrix(0.5*DRotAngle)
            );
        }
    }
}

tmp<surfaceVectorField> rotationAngleDerivative
(
    const volTensorField& rotMat
)
{
    tmp<surfaceVectorField> tRotAngDeriv
    (    
        new surfaceVectorField
        (
            IOobject
            (
                "rotAngleDerivative",
                rotMat.time().timeName(),
                rotMat.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rotMat.mesh(),
            dimensionedVector("zero", dimless/dimLength, vector::zero)
        )
    );

    vectorField& snGradRotAngleI = tRotAngDeriv().internalField();

    const tensorField& RI = rotMat.internalField();
    
    const fvMesh& mesh = rotMat.mesh();
    const unallocLabelList& own = mesh.owner();
    const unallocLabelList& nei = mesh.neighbour();
    const scalarField& deltaCoeffsI = mesh.deltaCoeffs();

    forAll(snGradRotAngleI, faceI)
    {
        tensor DR = (RI[own[faceI]].T() & RI[nei[faceI]]);

        vector DRotAngle = pseudoVector(DR);

        snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];
    }

    forAll(tRotAngDeriv().boundaryField(), patchI)
    {
        const tensorField& RI = rotMat.boundaryField()[patchI];

        vectorField& snGradRotAngleI =
            tRotAngDeriv().boundaryField()[patchI];

        const labelList& own =
            mesh.boundary()[patchI].faceCells();
        const scalarField& deltaCoeffsI =
            mesh.deltaCoeffs().boundaryField()[patchI];

        forAll(snGradRotAngleI, faceI)
        {
            tensor DR = (rotMat[own[faceI]].T() & RI[faceI]);

            vector DRotAngle = pseudoVector(DR);

            snGradRotAngleI[faceI] = DRotAngle*deltaCoeffsI[faceI];
        }
    }
    
    return tRotAngDeriv;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
