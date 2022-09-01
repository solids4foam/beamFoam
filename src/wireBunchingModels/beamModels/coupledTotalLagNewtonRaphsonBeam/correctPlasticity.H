if (plasticity_)
{
    scalarField& CPI = CP_.internalField();

    scalarField& ePI = eP_.internalField();
    const scalarField& oldEPI = eP_.oldTime().internalField();

    surfaceVectorField Q = (CQ_ & (Gamma_ - GammaP_));
    const vectorField& QI = Q.internalField();

    surfaceVectorField M = (CM_ & (K_ - KP_));
    const vectorField& MI = M.internalField();

    vectorField& GammaPI = GammaP_.internalField();
    vectorField& KPI = KP_.internalField();
    
    const vectorField& oldGammaPI = GammaP_.oldTime().internalField();
    const vectorField& oldKPI = KP_.oldTime().internalField();

    diagTensor invD
    (
        1.0/EA_.value(),
        1.0/GA_.value(),
        1.0/EI_.value()
    );

    forAll(CPI, faceI)
    {
        scalar N = QI[faceI].x();
        // scalar Q = QI[faceI].y();
        scalar M = MI[faceI].z();
        scalar NP = A_.value()*(sigmaY0_.value() + EP_.value()*mag(ePI[faceI]));

        scalar Ftrial = sqr(N) + CNM_.value()*NP*mag(M) - sqr(NP);
        
        if (Ftrial > 0)
        {
            vector dEpsilonP
            (
                GammaPI[faceI].x() - oldGammaPI[faceI].x(),
                GammaPI[faceI].y() - oldGammaPI[faceI].y(),
                KPI[faceI].z() - oldKPI[faceI].z()
            );
            scalar dEP(ePI[faceI] - oldEPI[faceI]);

            vector DFDS(2*N, 0, CNM_.value()*sign(M)*NP);
            scalar DFDNP = CNM_.value()*mag(M) - 2*NP;

            vector RS = CPI[faceI]*DFDS - dEpsilonP;
            scalar RNP = CPI[faceI]*DFDNP - dEP;

            label m = 1;
            if (dEP < -SMALL)
            {
                m = -1;
            }

            label signM = 1;
            if (M < -SMALL)
            {
                signM = -1;
            }

            diagTensor M00
            (
                EA_.value()/(1.0 + 2*CPI[faceI]*EA_.value()),
                GA_.value(),
                EI_.value()
               *(1.0 + 2*CPI[faceI]*EP_.value()*A_.value()/m)
               /(
                    1.0
                  + (2*CPI[faceI] + sqr(CNM_.value())*EI_.value()*sqr(CPI[faceI]))
                   *EP_.value()*A_.value()/m
                )
            );

            vector M01
            (
                0,
                0,
                CPI[faceI]
               *(
                   signM*CPI[faceI]*EI_.value()*EP_.value()*A_.value()/m
                )
               /(
                    1.0
                  + (2*CPI[faceI] + sqr(CNM_.value())*EI_.value()*sqr(CPI[faceI]))
                   *EP_.value()*A_.value()/m
                )
            );

            scalar M11 =
               -(EP_.value()*A_.value()/m)
               /(
                    1.0
                  + (2*CPI[faceI] + sqr(CNM_.value())*EI_.value()*sqr(CPI[faceI]))
                   *EP_.value()*A_.value()/m
                );

            scalar B =
                (DFDS & ((M00 & RS) + (M01*RNP)))
              + (DFDNP*((M01 & RS) + (M11*RNP)));
            
            scalar C =
                (DFDS & ((M00 & DFDS) + (M01*DFDNP)))
              + (DFDNP*((M01 & DFDS) + (M11*DFDNP)));

            scalar DCP = (Ftrial - B)/C;

            vector DS =
            (
                (M00 & (-RS - DCP*DFDS))
              + (M01 * (-RNP - DCP*DFDNP))
            );
            
            scalar DNP =
            (
                (M01 & (-RS - DCP*DFDS))
              + (M11 * (-RNP - DCP*DFDNP))
            );

            vector DEpsilonP = -(invD & DS);
            scalar DEP = (m/(EP_.value()+SMALL))*DNP;

            CP_.value() += DCP;
            ePI[faceI] += DEP;
            GammaPI[faceI].x() += DEpsilonP.x();
            KPI[faceI].z() += DEpsilonP.z();
        }
    }
}
