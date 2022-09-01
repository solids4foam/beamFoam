if (true)
{
    Info << "Updating coefficients" << endl;

    surfaceVectorField DThetaf = fvc::interpolate(DTheta);  

    surfaceScalarField magDThetaf = mag(DThetaf) + SMALL;

    surfaceTensorField DThetaHat = spinTensor(DThetaf);

    DLambda = II + (Foam::sin(magDThetaf)/magDThetaf)*DThetaHat
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*(DThetaHat & DThetaHat);
      // + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*(DThetaf*DThetaf)
      // - 2*sqr(Foam::sin(magDThetaf/2))*II;

    surfaceTensorField DT =
        (Foam::sin(magDThetaf)/magDThetaf)*II
      + ((1.0-Foam::sin(magDThetaf)/magDThetaf)/sqr(magDThetaf))
       *(DThetaf*DThetaf)
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*DThetaHat;

    CQDW = (DLambda & (CQ & DLambda.T()));
    
    surfaceTensorField A =
        (Foam::sin(magDThetaf)/magDThetaf)*II
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*DThetaHat; 
    
    CQDTheta = -(CQDW & (A & mulByPermutationTensor(dRuDs)))
      + (A & mulByPermutationTensor(Q.oldTime()));
    
    CMDTheta = ((DLambda & (CM & DLambda.T())) & DT);
    // CMDTheta = (DLambda & (CM & DT));
    CMDTheta2 = (A & mulByPermutationTensor(M.oldTime()));

    dQ.storePrevIter();
    // dQ = -(CQDW & (DThetaHat & dRuDs)) - (CQDTheta & DThetaf);
    dQ = (CQDW & ((II - DLambda) & dRuDs))
      + ((DLambda - II) & Q.oldTime())
      - (CQDTheta & DThetaf);
    dQ.relax();

    dM.storePrevIter();
    dM = ((DLambda - II) & M.oldTime())
      - (CMDTheta2 & DThetaf);
    dM.relax();
}
