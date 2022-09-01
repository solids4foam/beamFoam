{
    #include "updateRotation.H"

    const surfaceVectorField DWf = fvc::interpolate(DW);
    const vectorField& DWfI = DWf.internalField();

    const tensorField& DLambdaI = DLambda.internalField();
    const faceList& faces = mesh.faces();

    const vectorField& points = mesh.points();

    vectorField& pointDWI = pointDW.internalField();
    
    forAll(DWfI, faceI)
    {
        const face& curFace = faces[faceI];

        vector C0 = curFace.centre(points);

        forAll(curFace, pointI)
        {
            label curPoint = curFace[pointI];

            vector oldR = points[curPoint] - C0;

            vector newR = C0 + DWfI[faceI] + (DLambdaI[faceI] & oldR);

            pointDWI[curPoint] = newR - points[curPoint];
        }
    }

    forAll(DWf.boundaryField(), patchI)
    {
        const vectorField& pDWf =
            DWf.boundaryField()[patchI];

        const tensorField& pDLambda =
            DLambda.boundaryField()[patchI];
        
        const label start =
            mesh.boundaryMesh()[patchI].start();

        forAll(pDWf, faceI)
        {
            const face& curFace = faces[start + faceI];
            vector C0 = curFace.centre(points);

            // Info << C0 << endl;
            
            forAll(curFace, pointI)
            {
                label curPoint = curFace[pointI];

                vector oldR = points[curPoint] - C0;

                vector newR = C0 + pDWf[faceI] + (pDLambda[faceI] & oldR);

                pointDWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    vectorField newPoints = points + pointDWI;

    mesh.movePoints(newPoints);

    Ru = mesh.C();
    dRuDs = fvc::snGrad(Ru);
    dRuDs /= mag(dRuDs);

#   include "correctCurvedBeam.H"

    // dRuDs.boundaryField()[0] = vector(-1, 0, 0);

    // DTheta = dimensionedVector("0", DTheta.dimensions(), vector::zero);
    // Info << DTheta << endl;
}
