{
    const vectorField& refWfI = refWf.primitiveField();
    const tensorField& refLambdafI = refLambdaf.primitiveField();
    const faceList& faces = mesh.faces();

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.constant(),
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

    vectorField& pointWI = pointW.primitiveFieldRef();

    forAll(refWfI, faceI)
    {
        const face& curFace = faces[faceI];
        const vector C0 = curFace.centre(points);

        forAll(curFace, pointI)
        {
            const label curPoint = curFace[pointI];
            const vector oldR = points[curPoint] - C0;

            const vector newR =
                C0 + refWfI[faceI] + (refLambdafI[faceI] & oldR);

            pointWI[curPoint] = newR - points[curPoint];
        }
    }

    forAll(refWf.boundaryField(), patchI)
    {
        const vectorField& pRefWf =
            refWf.boundaryField()[patchI];

        const tensorField& pRefLambdaf =
            refLambdaf.boundaryField()[patchI];

        const label start =
            mesh.boundaryMesh()[patchI].start();

        forAll(pRefWf, faceI)
        {
            const face& curFace = faces[start + faceI];
            const vector C0 = curFace.centre(points);

            forAll(curFace, pointI)
            {
                const label curPoint = curFace[pointI];
                const vector oldR = points[curPoint] - C0;

                const vector newR =
                    C0 + pRefWf[faceI] + (pRefLambdaf[faceI] & oldR);

                pointWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    Info<< "Writing reference point displacement field" << endl;
    pointW.write();
    Info<< "Done setting circular arc reference point displacement" << endl;
}
