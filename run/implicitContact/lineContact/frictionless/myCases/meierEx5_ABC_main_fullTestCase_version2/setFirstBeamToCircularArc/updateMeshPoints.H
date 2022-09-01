{
    const vectorField& refWfI = refWf.internalField();

    const tensorField& refLambdaI = refLambda.internalField();
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
	
    pointMesh pMesh(mesh);

    pointVectorField pointW
    (
        IOobject
        (
            "pointW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector("0", dimLength, vector::zero)
    );

    vectorField& pointWI = pointW.internalField();


   // vectorField pointDWI(mesh.nPoints(), vector::zero);
    
    forAll(refWfI, faceI)
    {
        const face& curFace = faces[faceI];

        vector C0 = curFace.centre(points);

        forAll(curFace, pointI)
        {
            label curPoint = curFace[pointI];

            vector oldR = points[curPoint] - C0;

            vector newR = C0 + refWfI[faceI] + (refLambdaI[faceI] & oldR);

            pointWI[curPoint] = newR - points[curPoint];
        }
    }

    forAll(refWf.boundaryField(), patchI)
    {
        const vectorField& pRefWf =
            refWf.boundaryField()[patchI];

        const tensorField& pRefLambda =
            refLambda.boundaryField()[patchI];
        
        const label start =
            mesh.boundaryMesh()[patchI].start();

        forAll(pRefWf, faceI)
        {
            const face& curFace = faces[start + faceI];
            vector C0 = curFace.centre(points);

            forAll(curFace, pointI)
            {
                label curPoint = curFace[pointI];

                vector oldR = points[curPoint] - C0;

                vector newR = C0 + pRefWf[faceI] + (pRefLambda[faceI] & oldR);

                pointWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    // points = points + pointDWI;

    // Info << "Writing deformed mesh points ...";
    // points.write();
    // Info << "done" << endl;
    
    Info << "Writing reference point displacement field";
    pointW.write();
    Info << "done" << endl;
}
