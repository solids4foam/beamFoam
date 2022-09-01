if (true)
{
    // Info << "Updating coefficients" << endl;

    surfaceVectorField DThetaf = fvc::interpolate(DTheta_);  

    surfaceScalarField magDThetaf = mag(DThetaf) + SMALL;

    surfaceTensorField DThetaHat = spinTensor(DThetaf);

    dimensionedTensor II("I", dimless, tensor::I);

    DLambda_.storePrevIter();
    DLambda_ = II + (Foam::sin(magDThetaf)/magDThetaf)*DThetaHat
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*(DThetaHat & DThetaHat);
      // + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*(DThetaf*DThetaf)
      // - 2*sqr(Foam::sin(magDThetaf/2))*II;

    // scalar tmpDTheta = DTheta_.boundaryField()[0][0].z();
    
    // tensor tmpDLambda
    //     (
    //      ::cos(tmpDTheta), -::sin(tmpDTheta), 0,
    //      ::sin(tmpDTheta),  ::cos(tmpDTheta), 0,
    //         0, 0, 1
    //     );

    // DLambda_ = dimensionedTensor("DLambda", dimless, tmpDLambda);
    
    DLambda_.relax();

    Lambda_ = (DLambda_ & Lambda_.oldTime());
    
    
    // Info << average(DLambda_) << endl;
    
    surfaceTensorField DT =
        (Foam::sin(magDThetaf)/magDThetaf)*II
      + ((1.0-Foam::sin(magDThetaf)/magDThetaf)/sqr(magDThetaf))
       *(DThetaf*DThetaf)
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*DThetaHat;

    // Info << average(DT) << endl;
    
    CQDW_ = (Lambda_ & (CQ_ & Lambda_.T()));
    // CQDW_ = (DLambda_ & (CQ_ & DLambda_.T()));

    surfaceTensorField A =
        (Foam::sin(magDThetaf)/magDThetaf)*II
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*DThetaHat; 

    CQDTheta_ = -(CQDW_ & (A & mulByPermutationTensor(dRuDs_)))
      + (A & mulByPermutationTensor(Q_.oldTime()));

    // CMDTheta_ = (DLambda_ & (CM_ & DT));
    CMDTheta_ = ((Lambda_ & (CM_ & Lambda_.T())) & DT);

    CMDTheta2_ = (A & mulByPermutationTensor(M_.oldTime()));
}
