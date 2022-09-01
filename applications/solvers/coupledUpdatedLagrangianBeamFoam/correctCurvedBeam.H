// Correct cell centres for curved beam
{
    const surfaceVectorField& Cf = mesh.Cf();
    const vectorField& CfI = Cf.internalField();

    label nPoints = mesh.nCells()+1;

    vectorField beamPoints(nPoints, vector::zero);
    beamPoints[0] = Cf.boundaryField()[startPatchIndex][0];
    beamPoints[nPoints-1] = Cf.boundaryField()[endPatchIndex][0];
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
    // vectorField oldCentres = mesh.cellCentres();
    const_cast<vectorField&>(mesh.cellCentres()) = newCentres;

    // Info << "delta, max: " << max(mag(mesh.cellCentres() - oldCentres))
    //      << ", avg: " << average(mag(mesh.cellCentres() - oldCentres))
    //      << ", min: " << min(mag(mesh.cellCentres() - oldCentres)) << endl;

    Ru = mesh.C();
    dRuDs = fvc::snGrad(Ru);
    dRuDs /= mag(dRuDs);

    // vectorField dRdS = spline.dRdS();
    // for (label i=0; i<dRuDs.size(); i++)
    // {  
    //     dRuDs[i] = dRdS[i+1];
    // }
    // dRuDs.boundaryField()[startPatchIndex][0] = -dRdS[0];
    // dRuDs.boundaryField()[endPatchIndex][0] = dRdS[dRdS.size()-1];
}
