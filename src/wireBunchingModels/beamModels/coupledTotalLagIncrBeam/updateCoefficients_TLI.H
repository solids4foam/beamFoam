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
    DLambda_.relax();

    DT_ = (Foam::sin(magDThetaf)/magDThetaf)*II
      + ((1.0-Foam::sin(magDThetaf)/magDThetaf)/sqr(magDThetaf))
       *(DThetaf*DThetaf)
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*DThetaHat;

    Lambda_ = (DLambda_ & Lambda_.oldTime());
    CQDW_ = (Lambda_ & (CQ_ & Lambda_.T()));

    surfaceTensorField A =
        (Foam::sin(magDThetaf)/magDThetaf)*II
      + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*DThetaHat;

    CQDTheta_ =
       -(
            CQDW_
          & (
                A
              & mulByPermutationTensor
                (
                    (refLambda_ & i_)
                  + fvc::snGrad(W_.oldTime())
                )
            )
        )
      + (
            A
          & mulByPermutationTensor
            (
                Lambda_.oldTime()
              & (
                    CQ_
                  & (
                        gamma_.oldTime()
                    )
                )
            )
        );
    // CQTheta_ = -(CQW_ & (A & mulByPermutationTensor(dRuDs_)));

    CMDTheta_ = (Lambda_ & (CM_ & (Lambda_.oldTime().T() & DT_.T())));
    // CMTheta_ = ((Lambda_ & (CM_ & Lambda_.T())) & T);

    CMDTheta2_ =
    (
        A
      & mulByPermutationTensor
        (
            Lambda_.oldTime()
          & (
                CM_
              & (
                    kappa_.oldTime()
                )
            )
        )
    );

    // Explicit Q and M
    {
        // surfaceVectorField refThetaf = fvc::interpolate(refTheta_);
        // surfaceScalarField magRefThetaf = mag(refThetaf) + SMALL;
        // surfaceTensorField refThetaHat = spinTensor(refThetaf);

        // refLambda_ = II + (Foam::sin(magRefThetaf)/magRefThetaf)*refThetaHat
        //   + ((1.0-Foam::cos(magRefThetaf))/sqr(magRefThetaf))
        //    *(refThetaHat & refThetaHat);

        // surfaceTensorField refT =
        //     (Foam::sin(magRefThetaf)/magRefThetaf)*II
        //   + ((1.0-Foam::sin(magRefThetaf)/magRefThetaf)/sqr(magRefThetaf))
        //    *(refThetaf*refThetaf)
        //   + ((1.0-Foam::cos(magRefThetaf))/sqr(magRefThetaf))*refThetaHat;

        // explicitQ_ =
        // (
        //     Q_.oldTime()
        // );

        explicitQ_ =
        (
            Lambda_.oldTime()
          & (
                CQ_
              & (
                    gamma_.oldTime()
                )
            )
        );

        // explicitQ_ =
        // (
        //     (
        //         Lambda_
        //       & (
        //             CQ_
        //           & (
        //                 gamma_.oldTime()
        //             )
        //         )
        //     )
        // );

        // explicitQ_ = (CQW_ & ((refLambda_ - II) &  i_));

        // explicitM_ =
        // (
        //     M_.oldTime()
        // );

        explicitM_ =
        (
            (
                Lambda_.oldTime()
              & (
                    CM_
                  & (
                        kappa_.oldTime()
                    )
                )
            )
          + (
                Lambda_
              & (
                    CM_
                  & (
                      - (refT_.T() & fvc::snGrad(refTheta_))
                    )
                )
            )
        );

        // explicitM_ = -((Lambda_ & (CM_ & refT)) & fvc::snGrad(refTheta_));
    }
}
