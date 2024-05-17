surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

if (true)
{
  // Info << "Updating coefficients" << endl;

  // Total rotation matrix
  surfaceTensorField Lambdaf((Lambdaf_ & refLambdaf_));

  CQW_ = (Lambdaf & (CQ_ & Lambdaf.T()));

  explicitQ_ =  (Lambdaf & (CQ_ & (Gamma_)));

  explicitM_ =  (Lambdaf & (CM_ & (K_)));

  CQTheta_ =
  (
      (
          Lambdaf & (CQ_ & Lambdaf.T() )
      )
      & spinTensor(dRdS)
  )
  - spinTensor(Q_);
  // - spinTensor(explicitQ_);

  CQDTheta_ = (Lambdaf & (CDQDK_ & Lambdaf.T())); // Check for Kirchhoff beam

  CMTheta_ = ( Lambdaf & (CM_ & Lambdaf.T()) );

    // CMTheta2_ = -spinTensor(Lambdaf & (CM_ & (K_ - KP_)))

   CMTheta2_ = -spinTensor(M_);
   //CMTheta2_ = -spinTensor(explicitM_);

    CMQW_ =
        0.5
       *(
            (
                spinTensor(dRdS)
              & ( Lambdaf & (CQ_ & Lambdaf.T()) )
            )
          - spinTensor(Q_)
          // - spinTensor(explicitQ_)
        )/mesh().deltaCoeffs();
      // + (Lambdaf & (CDMDGamma_ & Lambdaf.T())); // Check for Kirchhoff beam

    CMQTheta_ =
        0.5
       *(
            (
                (
                    spinTensor(dRdS)
                  & ( Lambdaf & (CQ_ & Lambdaf.T()) )
                )
              & spinTensor(dRdS)
            )
          - (
                spinTensor(dRdS)
              & (
                    spinTensor(Q_)
                  //  spinTensor(explicitQ_)
                )
            )
        )
       /mesh().deltaCoeffs();

    explicitMQ_ =
        0.5*(spinTensor(dRdS) & explicitQ_)
       /mesh().deltaCoeffs();

    // Correct at boundary
    forAll(CMQW_.boundaryFieldRef(), patchI)
    {
        if (!CMQW_.boundaryFieldRef()[patchI].coupled())
        {
            CMQW_.boundaryFieldRef()[patchI] *= 2;
            CMQTheta_.boundaryFieldRef()[patchI] *= 2;
            explicitMQ_.boundaryFieldRef()[patchI] *= 2;
        }
    }
}
