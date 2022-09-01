    const dictionary& beamControl =
        beam().mesh().solutionDict().subDict("beamFoam");

    Foam::interpolationTable<Foam::scalar> deltaTseries;
    if (beamControl.found("deltaTseries"))
    {
        Info<< "deltaT is time-varying" << endl;
        deltaTseries =
            interpolationTable<scalar>
            (
                beamControl.subDict("deltaTseries")
            );
    }


    // int nCorr(readInt(stressControl.lookup("nCorrectors")));
    // scalar convergenceTolerance(readScalar(stressControl.lookup("U")));
