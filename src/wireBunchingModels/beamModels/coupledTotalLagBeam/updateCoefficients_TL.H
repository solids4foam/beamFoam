if (true)
{
    // Info << "Updating coefficients" << endl;

    surfaceVectorField Thetaf = fvc::interpolate(Theta_);  

    surfaceScalarField magThetaf = mag(Thetaf) + SMALL;

    surfaceTensorField ThetaHat = spinTensor(Thetaf);

    dimensionedTensor II("I", dimless, tensor::I);

    Lambda_.storePrevIter();
    Lambda_ = II + (Foam::sin(magThetaf)/magThetaf)*ThetaHat
      + ((1.0-Foam::cos(magThetaf))/sqr(magThetaf))*(ThetaHat & ThetaHat);
    Lambda_.relax();

    // Check pseudo vector calculation
    if (false)
    {
        vector theta = Theta_.boundaryField()[endPatchIndex()][0];
        tensor rotationMatrix = Lambda_.boundaryField()[endPatchIndex()][0];
        vector thetaPrime = pseudoVector(rotationMatrix);

        Info << theta << ", " << thetaPrime << endl;
    }

    surfaceTensorField T =
        (Foam::sin(magThetaf)/magThetaf)*II
      + ((1.0-Foam::sin(magThetaf)/magThetaf)/sqr(magThetaf))*(Thetaf*Thetaf)
      + ((1.0-Foam::cos(magThetaf))/sqr(magThetaf))*ThetaHat;

    CQW_ = (Lambda_ & (CQ_ & Lambda_.T()));

    surfaceTensorField A =
        (Foam::sin(magThetaf)/magThetaf)*II
      + ((1.0-Foam::cos(magThetaf))/sqr(magThetaf))*ThetaHat;

    CQTheta_ = -(CQW_ & (A & mulByPermutationTensor(i_)));
    // CQTheta_ = -(CQW_ & (A & mulByPermutationTensor(dRuDs_)));

    CMTheta_ = (Lambda_ & (CM_ & T.T()));
    // CMTheta_ = ((Lambda_ & (CM_ & Lambda_.T())) & T);

    // Explicit Q and M
    {
        surfaceVectorField refThetaf = fvc::interpolate(refTheta_);  

        surfaceScalarField magRefThetaf = mag(refThetaf) + SMALL;

        surfaceTensorField refThetaHat = spinTensor(refThetaf);

        refLambda_ = II + (Foam::sin(magRefThetaf)/magRefThetaf)*refThetaHat
          + ((1.0-Foam::cos(magRefThetaf))/sqr(magRefThetaf))
           *(refThetaHat & refThetaHat);

        surfaceTensorField refT =
            (Foam::sin(magRefThetaf)/magRefThetaf)*II
          + ((1.0-Foam::sin(magRefThetaf)/magRefThetaf)/sqr(magRefThetaf))
           *(refThetaf*refThetaf)
          + ((1.0-Foam::cos(magRefThetaf))/sqr(magRefThetaf))*refThetaHat;

        explicitQ_ = -(CQW_ & ((II - refLambda_) &  i_));
        explicitM_ = -((Lambda_ & (CM_ & refT.T())) & fvc::snGrad(refTheta_));
    }
}
