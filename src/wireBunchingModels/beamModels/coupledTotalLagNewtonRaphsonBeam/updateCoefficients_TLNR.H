surfaceVectorField dRdS = dR0Ds_ + fvc::snGrad(W_);

if (true)
{
    // Info << "Updating coefficients" << endl;

    // Total rotation matrix
    surfaceTensorField Lambda = (Lambda_ & refLambda_);

    CQW_ = (Lambda & (CDQ_ & Lambda.T()));
    // CQW_ = (Lambda & (CDQ_ & (kirchhoffTransTensor_ & Lambda.T())));

    // Pout << Lambda.boundaryField() << endl;

    // W_.correctBoundaryConditions();

    // // Update cell-centre mean line tangent
    // if (true)
    // {
    //     volVectorField curCentres = mesh().C() + W_;
    //     dimensionedVector i("one", dimless, vector(1, 0, 0));
    //     dRdS_ = (i & fvc::grad(curCentres));
    //     dRdS_ /= mag(dRdS_) + SMALL;
    // }

    //Info << Gamma_ << endl;
	
    explicitQ_ =  (Lambda & (CQ_ & (Gamma_ - GammaP_)));
    explicitM_ =  (Lambda & (CM_ & (K_ - KP_)));

    if (updatedLagrangian_) // Check for Kirchhoff beam and plasticity
    {
        explicitQ_ += (Lambda_ & Q_.oldTime());
        explicitM_ += (Lambda_ & M_.oldTime());
    }

    // #include "correctBoundaryForcesAndMoments.H"

    CQTheta_ =
        (
            (
                Lambda & (CDQ_ & Lambda.T() )
            )
          & spinTensor(dRdS)
        )
       - spinTensor(Q_);
     // - spinTensor(explicitQ_);

    CQDTheta_ = (Lambda & (CDQDK_ & Lambda.T())); // Check for Kirchhoff beam

    CMTheta_ = ( Lambda & (CDM_ & Lambda.T()) );

    // CMTheta2_ = -spinTensor(Lambda & (CM_ & (K_ - KP_)))
    
   CMTheta2_ = -spinTensor(M_);
   //CMTheta2_ = -spinTensor(explicitM_);
    
      // + ((Lambda & (CDMDGamma_ & Lambda.T())) & spinTensor(dRdS));

    CMQW_ =
        0.5
       *(
            (
                spinTensor(dRdS)
              & ( Lambda & (CDQ_ & Lambda.T()) )
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
                  & ( Lambda & (CDQ_ & Lambda.T()) )
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

    // if (updatedLagrangian_) // Check for Kirchhoff beam and plasticity
    // {
        // CQTheta_ -= spinTensor(Lambda_ & Q_.oldTime());
        // CMTheta2_ -= spinTensor(Lambda_ & M_.oldTime());
    
        // CMQW_ -= 0.5*spinTensor(Lambda_ & Q_.oldTime())
        //    /mesh().deltaCoeffs();
    
        // CMQTheta_ -=
        //     0.5*(spinTensor(dRdS) & spinTensor(Lambda_ & Q_.oldTime()))
        //    /mesh().deltaCoeffs();
    
        // explicitQ_ += (Lambda_ & Q_.oldTime());
        // explicitM_ += (Lambda_ & M_.oldTime());
        // explicitMQ_ += 0.5*(spinTensor(dRdS) & (Lambda_ & Q_.oldTime()))
        //    /mesh().deltaCoeffs();
    // }

    // Correct at boundary
    forAll(CMQW_.boundaryField(), patchI)
    {
        if (!CMQW_.boundaryField()[patchI].coupled())
        {
            CMQW_.boundaryField()[patchI] *= 2;
            CMQTheta_.boundaryField()[patchI] *= 2;
            explicitMQ_.boundaryField()[patchI] *= 2;
        }
    }

    // // Correct drds x Q term
    // explicitMQ_ = (Lambda.T() & explicitMQ_);
    // explicitMQ_.replace(0, 0);
    // explicitMQ_ = (Lambda & explicitMQ_);

    
    // Info << "Updating coefficients" << endl;
}
