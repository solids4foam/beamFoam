{
    const vectorField& DWfI = DWf.internalField();

    const tensorField& DLambdaI = DLambda.internalField();
    const faceList& faces = mesh.faces();

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.constant(),
            // runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    vectorField pointDWI(mesh.nPoints(), vector::zero);
    
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

            forAll(curFace, pointI)
            {
                label curPoint = curFace[pointI];

                vector oldR = points[curPoint] - C0;

                vector newR = C0 + pDWf[faceI] + (pDLambda[faceI] & oldR);

                pointDWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    points = points + pointDWI;

    Info << "Writing deformed mesh points ...";
    points.write();
    Info << "done" << endl;
}
