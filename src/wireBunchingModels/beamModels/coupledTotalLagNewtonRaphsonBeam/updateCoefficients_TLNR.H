surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

if (true)
{
  // Info << "Updating coefficients" << endl;

  // Total rotation matrix
  surfaceTensorField Lambda((Lambda_ & refLambda_));

  CQW_ = (Lambda & (CQ_ & Lambda.T()));

  explicitQ_ =  (Lambda & (CQ_ & (Gamma_)));

  explicitM_ =  (Lambda & (CM_ & (K_)));

  CQTheta_ =
  (
      (
          Lambda & (CQ_ & Lambda.T() )
      )
      & spinTensor(dRdS)
  )
  - spinTensor(Q_);
  // - spinTensor(explicitQ_);

  CQDTheta_ = (Lambda & (CDQDK_ & Lambda.T())); // Check for Kirchhoff beam

  CMTheta_ = ( Lambda & (CM_ & Lambda.T()) );

    // CMTheta2_ = -spinTensor(Lambda & (CM_ & (K_ - KP_)))

   CMTheta2_ = -spinTensor(M_);
   //CMTheta2_ = -spinTensor(explicitM_);

    CMQW_ =
        0.5
       *(
            (
                spinTensor(dRdS)
              & ( Lambda & (CQ_ & Lambda.T()) )
            )
          - spinTensor(Q_)
          // - spinTensor(explicitQ_)
        )/mesh().deltaCoeffs();
      // + (Lambda & (CDMDGamma_ & Lambda.T())); // Check for Kirchhoff beam

    CMQTheta_ =
        0.5
       *(
            (
                (
                    spinTensor(dRdS)
                  & ( Lambda & (CQ_ & Lambda.T()) )
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
