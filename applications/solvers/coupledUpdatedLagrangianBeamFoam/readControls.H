const dictionary& solidControl =
    mesh.solutionDict().subDict("solidMechanics");

const int nCorr(solidControl.lookupOrDefault<int>("nCorrectors", 1000));

const scalar convergenceTol
(
    solidControl.lookupOrDefault<scalar>("convergenceTol", 1e-6)
);


