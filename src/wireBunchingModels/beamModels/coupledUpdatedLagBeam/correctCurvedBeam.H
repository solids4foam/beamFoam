// Correct cell centres for curved beam
if (true)
{
    dynamicFvMesh& mesh =
        const_cast<dynamicFvMesh&>(this->mesh());


    // mesh.movePoints(mesh.points());

    // vectorField newCentres = mesh.cellCentres();
    // newCentres.replace(1, 0);
    // newCentres.replace(2, 0);
    // const_cast<vectorField&>(mesh.cellCentres()) = newCentres;

    // vectorField newFaceCentres = mesh.faceCentres();
    // for (label i=0; i<mesh.nInternalFaces(); i++)
    // {
    //     newFaceCentres[i].y() = 0;
    //     newFaceCentres[i].z() = 0;
    // }
    
    // label firstPatchFace = mesh.boundaryMesh()[startPatchIndex()].start();
    // newFaceCentres[firstPatchFace].y() = 0;
    // newFaceCentres[firstPatchFace].z() = 0;

    // label lastPatchFace = mesh.boundaryMesh()[endPatchIndex()].start();
    // newFaceCentres[lastPatchFace].y() = 0;
    // newFaceCentres[lastPatchFace].z() = 0;
        
    // const_cast<vectorField&>(mesh.faceCentres()) = newFaceCentres;


    const surfaceVectorField& Cf = mesh.Cf();
    const vectorField& CfI = Cf.internalField();

    label nPoints = mesh.nCells()+1;

    vectorField beamPoints(nPoints, vector::zero);
    beamPoints[0] = Cf.boundaryField()[startPatchIndex()][0];
    beamPoints[nPoints-1] = Cf.boundaryField()[endPatchIndex()][0];
    for (label i=0; i<CfI.size(); i++)
    {
        beamPoints[i+1] = Cf[i];
    }

    cubicSpline spline
    (
        beamPoints,
        cubicSpline::CLAMPED_CALC,
        vector(0, 0, 0),
        cubicSpline::CLAMPED_CALC,
        vector(0, 0, 0)
    );

    mesh.movePoints(mesh.points());
    vectorField newCentres = spline.midPoints();
    const_cast<vectorField&>(mesh.cellCentres()) = newCentres;

    // Set segment lengths
    this->L().internalField() = spline.segLengths();
}
