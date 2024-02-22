    {
        const tensorField& pCQW = CQW_.boundaryField()[patchI];
        const tensorField& pCQTheta = CQTheta_.boundaryField()[patchI];
        const tensorField& pCQDTheta = CQDTheta_.boundaryField()[patchI];
        const tensorField& pCMTheta = CMTheta_.boundaryField()[patchI];
        const tensorField& pCMTheta2 = CMTheta2_.boundaryField()[patchI];
        const tensorField& pCMQW = CMQW_.boundaryField()[patchI];
        const tensorField& pCMQTheta = CMQTheta_.boundaryField()[patchI];

        const vectorField& pExplicitQ = explicitQ_.boundaryField()[patchI];
        const vectorField& pExplicitM = explicitM_.boundaryField()[patchI];
        const vectorField& pExplicitMQ = explicitMQ_.boundaryField()[patchI];

        const tensorField& pLambda = Lambda_.boundaryField()[patchI];

        ////// W equation

        if
        (
            isA<forceBeamDisplacementNRFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementNRFvPatchVectorField& pW =
                refCast<forceBeamDisplacementNRFvPatchVectorField>
                (
                    W_.boundaryFieldRef()[patchI]
                );
            // Source contribution
            forAll (pW, faceI)
            {
                source[fc[faceI]](0,0) -= pW.force()[faceI].x();
                source[fc[faceI]](1,0) -= pW.force()[faceI].y();
                source[fc[faceI]](2,0) -= pW.force()[faceI].z();

                // This part will be subtracted later
                source[fc[faceI]](0,0) += pExplicitQ[faceI].x();
                source[fc[faceI]](1,0) += pExplicitQ[faceI].y();
                source[fc[faceI]](2,0) += pExplicitQ[faceI].z();
            }

	    if
	    (
		isA<followerForceBeamDisplacementNRFvPatchVectorField>
		(
		    W_.boundaryField()[patchI]
		)
	    )
	    {
		vectorField followerForce = pW.force();

		tensorField ffDiag(-spinTensor(followerForce));


		forAll (pW, faceI)
		{
		    d[fc[faceI]](0,3) += ffDiag[faceI].xx();
		    d[fc[faceI]](0,4) += ffDiag[faceI].xy();
		    d[fc[faceI]](0,5) += ffDiag[faceI].xz();

		    d[fc[faceI]](1,3) += ffDiag[faceI].yx();
		    d[fc[faceI]](1,4) += ffDiag[faceI].yy();
		    d[fc[faceI]](1,5) += ffDiag[faceI].yz();

		    d[fc[faceI]](2,3) += ffDiag[faceI].zx();
		    d[fc[faceI]](2,4) += ffDiag[faceI].zy();
		    d[fc[faceI]](2,5) += ffDiag[faceI].zz();

		}
	    }
        }
        else if
        (
            isA<axialForceTransverseDisplacementNRFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementNRFvPatchVectorField& pW =
                refCast<axialForceTransverseDisplacementNRFvPatchVectorField>
                (
                    W_.boundaryFieldRef()[patchI]
                );

            scalarField pDelta
            (
                1.0/mesh().deltaCoeffs().boundaryField()[patchI]
            );
            vectorField tang
            (
            (
                Lambda_.boundaryField()[patchI]
              & dR0Ds_.boundaryField()[patchI]
            )
            );

            const vectorField& pWPrev
            (
                W_.prevIter().boundaryField()[patchI]
            );
            vectorField pWCorr(pW.refDisp() - pWPrev);

            tensorField pCQWt(pCQW - ((tang*tang) & pCQW));
            tensorField pCQThetat(pCQTheta - ((tang*tang) & pCQTheta));

            // Diag contribution
            forAll (pW, faceI)
            {
                d[fc[faceI]](0,0) += -pCQWt[faceI].xx()/pDelta[faceI];
                d[fc[faceI]](0,1) += -pCQWt[faceI].xy()/pDelta[faceI];
                d[fc[faceI]](0,2) += -pCQWt[faceI].xz()/pDelta[faceI];

                d[fc[faceI]](1,0) += -pCQWt[faceI].yx()/pDelta[faceI];
                d[fc[faceI]](1,1) += -pCQWt[faceI].yy()/pDelta[faceI];
                d[fc[faceI]](1,2) += -pCQWt[faceI].yz()/pDelta[faceI];

                d[fc[faceI]](2,0) += -pCQWt[faceI].zx()/pDelta[faceI];
                d[fc[faceI]](2,1) += -pCQWt[faceI].zy()/pDelta[faceI];
                d[fc[faceI]](2,2) += -pCQWt[faceI].zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pW, faceI)
            {
                vector WContrib =
                    (pCQWt[faceI] & pWCorr[faceI])/pDelta[faceI];

                source[fc[faceI]](0,0) -= WContrib.x();
                source[fc[faceI]](1,0) -= WContrib.y();
                source[fc[faceI]](2,0) -= WContrib.z();
            }

            if
            (
                isA<momentBeamRotationNRFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationNRFvPatchVectorField& pTheta
                (
                    refCast<momentBeamRotationNRFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    )
                );
                tensorField invCM(inv(pCMTheta/pDelta + pCMTheta2));

                // Diag contribution
                forAll (pTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQThetat[faceI]
                        & (invCM[faceI] & (pCMTheta[faceI]/pDelta[faceI]))
                    );

                    d[fc[faceI]](0,3) += CqCt.xx();
                    d[fc[faceI]](0,4) += CqCt.xy();
                    d[fc[faceI]](0,5) += CqCt.xz();

                    d[fc[faceI]](1,3) += CqCt.yx();
                    d[fc[faceI]](1,4) += CqCt.yy();
                    d[fc[faceI]](1,5) += CqCt.yz();

                    d[fc[faceI]](2,3) += CqCt.zx();
                    d[fc[faceI]](2,4) += CqCt.zy();
                    d[fc[faceI]](2,5) += CqCt.zz();
                }

                // Source contribution
                forAll (pTheta, faceI)
                {
                    tensor CqCt = pCQThetat[faceI];
                    vector thetaContrib =
                        (
                            CqCt
                          & (
                                invCM[faceI]
                              & (
                                    pTheta.moment()[faceI]
                                  - pExplicitM[faceI]
                                )
                            )
                        );

                    source[fc[faceI]](0,0) -= thetaContrib.x();
                    source[fc[faceI]](1,0) -= thetaContrib.y();
                    source[fc[faceI]](2,0) -= thetaContrib.z();
                }
            }
            else if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pTheta
                (
                    refCast<fixedValueFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    )
                );

                const vectorField& pThetaPrev =
                    Theta_.prevIter().boundaryField()[patchI];

                vectorField& pThetaCorr
                (
                    DTheta_.boundaryFieldRef()[patchI]
                );
                // pThetaCorrSta = (pLambda & pThetaCorrStar);
                // pThetaCorrStar = (pTheta - pThetaPrev);

                pThetaCorr = (pLambda & (pTheta - pThetaPrev));

                // Info << "(ax) pThetaCorr = " << pThetaCorr << endl;

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector thetaContrib
                    (
                        (pCQThetat[faceI] & pThetaCorr[faceI])
                    );
                    source[fc[faceI]](0,0) -= thetaContrib.x();
                    source[fc[faceI]](1,0) -= thetaContrib.y();
                    source[fc[faceI]](2,0) -= thetaContrib.z();
                }
            }

            // Add axial component of force
            vectorField pQa(tang*pW.axialForce());
            forAll (pQa, faceI)
            {
                source[fc[faceI]](0,0) -= pQa[faceI].x();
                source[fc[faceI]](1,0) -= pQa[faceI].y();
                source[fc[faceI]](2,0) -= pQa[faceI].z();
            }

            // Subtract explicit force contribution
            vectorField pExplicitQt(((tang*tang) & pExplicitQ));
            forAll (pExplicitQ, faceI)
            {
                source[fc[faceI]](0,0) += pExplicitQt[faceI].x();
                source[fc[faceI]](1,0) += pExplicitQt[faceI].y();
                source[fc[faceI]](2,0) += pExplicitQt[faceI].z();
            }
        }
        else if
        (
            isA<fixedValueFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pW
            (
                refCast<fixedValueFvPatchVectorField>
                (
                    W_.boundaryFieldRef()[patchI]
                )
            );

            const vectorField& pWPrev
            (
                W_.prevIter().boundaryField()[patchI]
            );
            vectorField& pWCorr(DW_.boundaryFieldRef()[patchI]);

            pWCorr = pW - pWPrev;

            scalarField pDelta
            (
                1.0/mesh().deltaCoeffs().boundaryField()[patchI]
            );
            // Diag contribution
            forAll (pW, faceI)
            {
                d[fc[faceI]](0,0) += -pCQW[faceI].xx()/pDelta[faceI];
                d[fc[faceI]](0,1) += -pCQW[faceI].xy()/pDelta[faceI];
                d[fc[faceI]](0,2) += -pCQW[faceI].xz()/pDelta[faceI];

                d[fc[faceI]](1,0) += -pCQW[faceI].yx()/pDelta[faceI];
                d[fc[faceI]](1,1) += -pCQW[faceI].yy()/pDelta[faceI];
                d[fc[faceI]](1,2) += -pCQW[faceI].yz()/pDelta[faceI];

                d[fc[faceI]](2,0) += -pCQW[faceI].zx()/pDelta[faceI];
                d[fc[faceI]](2,1) += -pCQW[faceI].zy()/pDelta[faceI];
                d[fc[faceI]](2,2) += -pCQW[faceI].zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pW, faceI)
            {
                vector WContrib =
                (
                    pCQW[faceI]
                  & (
                        pWCorr[faceI]
                      + vector(SMALL, SMALL, SMALL) // Not sure why
                    )
                )
               /pDelta[faceI];

                source[fc[faceI]](0,0) -= WContrib.x();
                source[fc[faceI]](1,0) -= WContrib.y();
                source[fc[faceI]](2,0) -= WContrib.z();
            }

            if
            (
                isA<momentBeamRotationNRFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationNRFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationNRFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                tensorField invCM(inv(pCMTheta/pDelta + pCMTheta2));

                // Diag contribution
                forAll (pTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQTheta[faceI]
                      & (invCM[faceI] & (pCMTheta[faceI]/pDelta[faceI]))
                    );

                    d[fc[faceI]](0,3) += CqCt.xx();
                    d[fc[faceI]](0,4) += CqCt.xy();
                    d[fc[faceI]](0,5) += CqCt.xz();

                    d[fc[faceI]](1,3) += CqCt.yx();
                    d[fc[faceI]](1,4) += CqCt.yy();
                    d[fc[faceI]](1,5) += CqCt.yz();

                    d[fc[faceI]](2,3) += CqCt.zx();
                    d[fc[faceI]](2,4) += CqCt.zy();
                    d[fc[faceI]](2,5) += CqCt.zz();
                }

                // Source contribution
                forAll (pTheta, faceI)
                {
                    tensor CqCt = pCQTheta[faceI];
                    vector thetaContrib =
                        (
                            CqCt
                          & (
                                invCM[faceI]
                              & (
                                    pTheta.moment()[faceI]
                                  - pExplicitM[faceI]
                                )
                            )
                        );

                    source[fc[faceI]](0,0) -= thetaContrib.x();
                    source[fc[faceI]](1,0) -= thetaContrib.y();
                    source[fc[faceI]](2,0) -= thetaContrib.z();
                }
            }
            else if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                // Info << Theta_.boundaryField()[patchI].type() << endl;

                const fixedValueFvPatchVectorField& pTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                const vectorField& pThetaPrev =
                    Theta_.prevIter().boundaryField()[patchI];

                vectorField& pThetaCorr =
                    DTheta_.boundaryFieldRef()[patchI];

                // pThetaCorrStar = (pTheta - pThetaPrev);
                pThetaCorr = (pLambda & (pTheta - pThetaPrev));

                // Info << "pThetaCorr = " << pThetaCorr << endl;

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector thetaContrib =
                    (
                        pCQTheta[faceI]
                      & pThetaCorr[faceI]
                    )
                  + (
                        pCQDTheta[faceI]
                      & pThetaCorr[faceI]
                    )
                   /pDelta[faceI];
                    // (pCQTheta[faceI] & pTheta[faceI]);

                    source[fc[faceI]](0,0) -= thetaContrib.x();
                    source[fc[faceI]](1,0) -= thetaContrib.y();
                    source[fc[faceI]](2,0) -= thetaContrib.z();
                }

                // Diag contribution
                forAll (pW, faceI)
                {
                    d[fc[faceI]](0,0) += -pCQDTheta[faceI].xx()/pDelta[faceI];
                    d[fc[faceI]](0,1) += -pCQDTheta[faceI].xy()/pDelta[faceI];
                    d[fc[faceI]](0,2) += -pCQDTheta[faceI].xz()/pDelta[faceI];

                    d[fc[faceI]](1,0) += -pCQDTheta[faceI].yx()/pDelta[faceI];
                    d[fc[faceI]](1,1) += -pCQDTheta[faceI].yy()/pDelta[faceI];
                    d[fc[faceI]](1,2) += -pCQDTheta[faceI].yz()/pDelta[faceI];

                    d[fc[faceI]](2,0) += -pCQDTheta[faceI].zx()/pDelta[faceI];
                    d[fc[faceI]](2,1) += -pCQDTheta[faceI].zy()/pDelta[faceI];
                    d[fc[faceI]](2,2) += -pCQDTheta[faceI].zz()/pDelta[faceI];
                }
            }
        }

        // Add explicit force conribution
        forAll (pExplicitQ, faceI)
        {
            source[fc[faceI]](0,0) -= pExplicitQ[faceI].x();
            source[fc[faceI]](1,0) -= pExplicitQ[faceI].y();
            source[fc[faceI]](2,0) -= pExplicitQ[faceI].z();
        }


        //////////// Theta equation

        if
        (
            isA<momentBeamRotationNRFvPatchVectorField>
            (
                Theta_.boundaryField()[patchI]
            )
        )
        {
            const momentBeamRotationNRFvPatchVectorField& pTheta =
                refCast<momentBeamRotationNRFvPatchVectorField>
                (
                    Theta_.boundaryFieldRef()[patchI]
                );

            // scalarField pDelta =
            //     1.0/mesh().deltaCoeffs().boundaryField()[patchI];

            // tensorField invCM = inv(pCMTheta/pDelta + pCMTheta2);

            // tensorField T0 = (I - (invCM & (pCMTheta/pDelta)));

            // tensorField pCMThetaPrime =
            //     ((pCMTheta/pDelta) & T0)
            //   - (pCMTheta2 & (invCM & (pCMTheta/pDelta)));

            // // Diag contribution from laplacian term
            // forAll (pTheta, faceI)
            // {
            //     d[fc[faceI]](3,3) += -pCMThetaPrime[faceI].xx();
            //     d[fc[faceI]](3,4) += -pCMThetaPrime[faceI].xy();
            //     d[fc[faceI]](3,5) += -pCMThetaPrime[faceI].xz();

            //     d[fc[faceI]](4,3) += -pCMThetaPrime[faceI].yx();
            //     d[fc[faceI]](4,4) += -pCMThetaPrime[faceI].yy();
            //     d[fc[faceI]](4,5) += -pCMThetaPrime[faceI].yz();

            //     d[fc[faceI]](5,3) += -pCMThetaPrime[faceI].zx();
            //     d[fc[faceI]](5,4) += -pCMThetaPrime[faceI].zy();
            //     d[fc[faceI]](5,5) += -pCMThetaPrime[faceI].zz();
            // }

            // Source contribution
            vectorField dM((pTheta.moment() - pExplicitM));
            forAll (pTheta, faceI)
            {
                source[fc[faceI]](3,0) -= dM[faceI].x();
                source[fc[faceI]](4,0) -= dM[faceI].y();
                source[fc[faceI]](5,0) -= dM[faceI].z();
            }
        }
        else if
        (
            isA<fixedValueFvPatchVectorField>
            (
                Theta_.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pTheta =
                refCast<fixedValueFvPatchVectorField>
                (
                    Theta_.boundaryFieldRef()[patchI]
                );

            const vectorField& pThetaPrev =
                Theta_.prevIter().boundaryField()[patchI];

            vectorField& pThetaCorr
            (
                DTheta_.boundaryFieldRef()[patchI]
            );
            // pThetaCorrStar = (pTheta - pThetaPrev);
            pThetaCorr = (pLambda & (pTheta - pThetaPrev));

            // Info << patchI << ", " << pThetaCorr << endl;

            scalarField pDelta
            (
                1.0/mesh().deltaCoeffs().boundaryField()[patchI]
            );
            // Diag contribution from laplacian term
            forAll (pTheta, faceI)
            {
                d[fc[faceI]](3,3) += -pCMTheta[faceI].xx()/pDelta[faceI];
                d[fc[faceI]](3,4) += -pCMTheta[faceI].xy()/pDelta[faceI];
                d[fc[faceI]](3,5) += -pCMTheta[faceI].xz()/pDelta[faceI];

                d[fc[faceI]](4,3) += -pCMTheta[faceI].yx()/pDelta[faceI];
                d[fc[faceI]](4,4) += -pCMTheta[faceI].yy()/pDelta[faceI];
                d[fc[faceI]](4,5) += -pCMTheta[faceI].yz()/pDelta[faceI];

                d[fc[faceI]](5,3) += -pCMTheta[faceI].zx()/pDelta[faceI];
                d[fc[faceI]](5,4) += -pCMTheta[faceI].zy()/pDelta[faceI];
                d[fc[faceI]](5,5) += -pCMTheta[faceI].zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pTheta, faceI)
            {
                // Source contribution from laplacian term
                source[fc[faceI]](3,0) -=
                    (
                        pCMTheta[faceI].xx()*pThetaCorr[faceI].x()
                      + pCMTheta[faceI].xy()*pThetaCorr[faceI].y()
                      + pCMTheta[faceI].xz()*pThetaCorr[faceI].z()
                    )
                   /pDelta[faceI];

                source[fc[faceI]](4,0) -=
                    (
                        pCMTheta[faceI].yx()*pThetaCorr[faceI].x()
                      + pCMTheta[faceI].yy()*pThetaCorr[faceI].y()
                      + pCMTheta[faceI].yz()*pThetaCorr[faceI].z()
                    )
                   /pDelta[faceI];

                source[fc[faceI]](5,0) -=
                    (
                        pCMTheta[faceI].zx()*pThetaCorr[faceI].x()
                      + pCMTheta[faceI].zy()*pThetaCorr[faceI].y()
                      + pCMTheta[faceI].zz()*pThetaCorr[faceI].z()
                    )
                   /pDelta[faceI];

                // Theta part
                source[fc[faceI]](3,0) -=
                    (
                        pCMTheta2[faceI].xx()*pThetaCorr[faceI].x()
                      + pCMTheta2[faceI].xy()*pThetaCorr[faceI].y()
                      + pCMTheta2[faceI].xz()*pThetaCorr[faceI].z()
                    );

                source[fc[faceI]](4,0) -=
                    (
                        pCMTheta2[faceI].yx()*pThetaCorr[faceI].x()
                      + pCMTheta2[faceI].yy()*pThetaCorr[faceI].y()
                      + pCMTheta2[faceI].yz()*pThetaCorr[faceI].z()
                    );

                source[fc[faceI]](5,0) -=
                    (
                        pCMTheta2[faceI].zx()*pThetaCorr[faceI].x()
                      + pCMTheta2[faceI].zy()*pThetaCorr[faceI].y()
                      + pCMTheta2[faceI].zz()*pThetaCorr[faceI].z()
                    );
            }
        }

        // Explicit moment
        forAll (pExplicitM, faceI)
        {
            source[fc[faceI]](3,0) -= pExplicitM[faceI].x();
            source[fc[faceI]](4,0) -= pExplicitM[faceI].y();
            source[fc[faceI]](5,0) -= pExplicitM[faceI].z();
        }


        // dr x Q term

        scalarField pDelta
        (
            1.0/mesh().deltaCoeffs().boundaryField()[patchI]
        );
        if
        (
            isA<forceBeamDisplacementNRFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementNRFvPatchVectorField& pW =
                refCast<forceBeamDisplacementNRFvPatchVectorField>
                (
                    W_.boundaryFieldRef()[patchI]
                );

            tensorField invCQW(inv(pCQW));
            // vectorField dwds0 = (invCQW & (pW.force() - pExplicitQ));

            vectorField pDRdS = dRdS.boundaryField()[patchI];

            // Source contribution
            forAll (pW, faceI)
            {
                // vector curSource = (pCMQW[faceI] & dwds0[faceI]);

                if (true) // This term is zero
                {
                    vector curSource =
                    (
                        spinTensor(pDRdS[faceI])
                      & (pW.force()[faceI] - pExplicitQ[faceI])
                    )
                   *pDelta[faceI];

                    source[fc[faceI]](3,0) -= curSource.x();
                    source[fc[faceI]](4,0) -= curSource.y();
                    source[fc[faceI]](5,0) -= curSource.z();
                }
            }

            if
            (
                isA<momentBeamRotationNRFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationNRFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationNRFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                if (false)
                {
                    tensorField A(pCQW);

                    tensorField B(pCMTheta);

                    tensorField C
                    (
                        B/pDelta
                      - spinTensor(pExplicitM)
                    );

                    tensorField D
                    (
                        (inv(A) & spinTensor(pExplicitQ))
                      - spinTensor(pDRdS)
                    );

                    vectorField src
                    (
                        (D & (inv(C) & (pTheta.moment() - pExplicitM)))
                      - (inv(A) & (pW.force() - pExplicitQ))
                    );
                    src = -(spinTensor(pW.force()) & src)*pDelta;

                    // src = -(spinTensor(pExplicitQ) & src)*pDelta;

                    tensorField diag((D & (inv(C) & B))/pDelta);
                    diag = -(spinTensor(pW.force()) & diag)*pDelta;
                    // diag = -(spinTensor(pExplicitQ) & diag)*pDelta;

                    // Source contribution
                    forAll (pTheta, faceI)
                    {
                        source[fc[faceI]](3,0) -= src[faceI].x();
                        source[fc[faceI]](4,0) -= src[faceI].y();
                        source[fc[faceI]](5,0) -= src[faceI].z();
                    }

                    // Diag contribution
                    forAll (diag, faceI)
                    {
                        d[fc[faceI]](3,3) += diag[faceI].xx();
                        d[fc[faceI]](3,4) += diag[faceI].xy();
                        d[fc[faceI]](3,5) += diag[faceI].xz();

                        d[fc[faceI]](4,3) += diag[faceI].yx();
                        d[fc[faceI]](4,4) += diag[faceI].yy();
                        d[fc[faceI]](4,5) += diag[faceI].yz();

                        d[fc[faceI]](5,3) += diag[faceI].zx();
                        d[fc[faceI]](5,4) += diag[faceI].zy();
                        d[fc[faceI]](5,5) += diag[faceI].zz();
                    }
                }

                if (true)
                {
                    scalarField pDelta
                    (
                        1.0/mesh().deltaCoeffs().boundaryField()[patchI]
                    );
                    tensorField invCM(inv(pCMTheta/pDelta + pCMTheta2));

                    vectorField A((invCM & (pTheta.moment() - pExplicitM)));
                    tensorField B((invCM & (pCMTheta/pDelta)));

                    vectorField src((pCMQTheta & A)
                        - (pCMQW & ((invCQW & pCQTheta) & A)));

                    tensorField diag((pCMQTheta & B)
                        - (pCMQW & ((invCQW & pCQTheta) & B)));

                    // Diag contribution
                    forAll (pTheta, faceI)
                    {
                        d[fc[faceI]](3,3) += diag[faceI].xx();
                        d[fc[faceI]](3,4) += diag[faceI].xy();
                        d[fc[faceI]](3,5) += diag[faceI].xz();

                        d[fc[faceI]](4,3) += diag[faceI].yx();
                        d[fc[faceI]](4,4) += diag[faceI].yy();
                        d[fc[faceI]](4,5) += diag[faceI].yz();

                        d[fc[faceI]](5,3) += diag[faceI].zx();
                        d[fc[faceI]](5,4) += diag[faceI].zy();
                        d[fc[faceI]](5,5) += diag[faceI].zz();
                    }

                    // Source contribution
                    forAll (pTheta, faceI)
                    {
                        source[fc[faceI]](3,0) -= src[faceI].x();
                        source[fc[faceI]](4,0) -= src[faceI].y();
                        source[fc[faceI]](5,0) -= src[faceI].z();
                    }
                }
            }

	    if
	    (
		isA<followerForceBeamDisplacementNRFvPatchVectorField>
		(
		    W_.boundaryField()[patchI]
		)
	    )
	    {
		scalarField pDelta
        (
                1.0/mesh().deltaCoeffs().boundaryField()[patchI]
        );
		vectorField pDRdS = dRdS.boundaryField()[patchI];

		vectorField followerForce = pW.force();

		//tensorField ffdRcrossQDiag(pW.size(),tensor::zero);


		tensorField ffdRcrossQDiag
        (
				-pDelta*
				(
				    spinTensor(pDRdS)
				    & spinTensor(followerForce)
				)
        );

		// Diag contribution
		forAll (pW, faceI)
		{
		    d[fc[faceI]](3,3) += ffdRcrossQDiag[faceI].xx();
		    d[fc[faceI]](3,4) += ffdRcrossQDiag[faceI].xy();
		    d[fc[faceI]](3,5) += ffdRcrossQDiag[faceI].xz();

		    d[fc[faceI]](4,3) += ffdRcrossQDiag[faceI].yx();
		    d[fc[faceI]](4,4) += ffdRcrossQDiag[faceI].yy();
		    d[fc[faceI]](4,5) += ffdRcrossQDiag[faceI].yz();

		    d[fc[faceI]](5,3) += ffdRcrossQDiag[faceI].zx();
		    d[fc[faceI]](5,4) += ffdRcrossQDiag[faceI].zy();
		    d[fc[faceI]](5,5) += ffdRcrossQDiag[faceI].zz();

		}
	    }

            else if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                const vectorField& pThetaPrev
                (
                    Theta_.prevIter().boundaryField()[patchI]
                );
                vectorField& pThetaCorr
                (
                    DTheta_.boundaryFieldRef()[patchI]
                );
                // pThetaCorrStar = (pTheta - pThetaPrev);
                pThetaCorr = (pLambda & (pTheta - pThetaPrev));

                vectorField dwds1(-((invCQW & pCQTheta) & pThetaCorr));

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector curSource =
                        (pCMQW[faceI] & dwds1[faceI])
                      + (pCMQTheta[faceI] & pThetaCorr[faceI]);

                    source[fc[faceI]](3,0) -= curSource.x();
                    source[fc[faceI]](4,0) -= curSource.y();
                    source[fc[faceI]](5,0) -= curSource.z();
                }
            }
        }
        else if
        (
            isA<axialForceTransverseDisplacementNRFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementNRFvPatchVectorField& pW =
                refCast<axialForceTransverseDisplacementNRFvPatchVectorField>
                (
                    W_.boundaryFieldRef()[patchI]
                );

            scalarField pDelta
            (
                1.0/mesh().deltaCoeffs().boundaryField()[patchI]
            );
            vectorField tang
            (
            (
                Lambda_.boundaryField()[patchI]
              & dR0Ds_.boundaryField()[patchI]
            )
            );

            // tensorField pCQWt = pCQW - ((tang*tang) & pCQW);
            tensorField pCQThetat(pCQTheta - ((tang*tang) & pCQTheta));

            tensorField invCQW(inv(pCQW));

            vectorField dwds0
            (
            (
                invCQW
              & (
                    pW.axialForce()*tang
                  - ((tang*tang) & pExplicitQ)
                )
            )
            );

            // Source contribution
            forAll (pW, faceI)
            {
                vector curSource = (pCMQW[faceI] & dwds0[faceI]);

                source[fc[faceI]](3,0) -= curSource.x();
                source[fc[faceI]](4,0) -= curSource.y();
                source[fc[faceI]](5,0) -= curSource.z();
            }

            if
            (
                isA<momentBeamRotationNRFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationNRFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationNRFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                scalarField pDelta
                (
                    1.0/mesh().deltaCoeffs().boundaryField()[patchI]
                );
                tensorField invCM(inv(pCMTheta/pDelta + pCMTheta2));

                vectorField A((invCM & (pTheta.moment() - pExplicitM)));
                tensorField B((invCM & (pCMTheta/pDelta)));

                vectorField src
                (
                 (pCMQTheta & A)
                 - (pCMQW & ((invCQW & pCQThetat) & A))
                );
                tensorField diag
                (
                    (pCMQTheta & B)
                  - (pCMQW & ((invCQW & pCQThetat) & B))
                );

                // Diag contribution
                forAll (pTheta, faceI)
                {
                    d[fc[faceI]](3,3) += diag[faceI].xx();
                    d[fc[faceI]](3,4) += diag[faceI].xy();
                    d[fc[faceI]](3,5) += diag[faceI].xz();

                    d[fc[faceI]](4,3) += diag[faceI].yx();
                    d[fc[faceI]](4,4) += diag[faceI].yy();
                    d[fc[faceI]](4,5) += diag[faceI].yz();

                    d[fc[faceI]](5,3) += diag[faceI].zx();
                    d[fc[faceI]](5,4) += diag[faceI].zy();
                    d[fc[faceI]](5,5) += diag[faceI].zz();
                }

                // Source contribution
                forAll (pTheta, faceI)
                {
                    source[fc[faceI]](3,0) -= src[faceI].x();
                    source[fc[faceI]](4,0) -= src[faceI].y();
                    source[fc[faceI]](5,0) -= src[faceI].z();
                }
            }
            else if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                const vectorField& pThetaPrev =
                    Theta_.prevIter().boundaryField()[patchI];

                vectorField& pThetaCorrStar
                (
                    DTheta_.boundaryFieldRef()[patchI]
                );
                pThetaCorrStar = (pTheta - pThetaPrev);
                vectorField pThetaCorr((pLambda & pThetaCorrStar));

                vectorField dwds1(-((invCQW & pCQThetat) & pThetaCorr));

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector curSource =
                        (pCMQW[faceI] & dwds1[faceI])
                      + (pCMQTheta[faceI] & pThetaCorr[faceI]);

                    source[fc[faceI]](3,0) -= curSource.x();
                    source[fc[faceI]](4,0) -= curSource.y();
                    source[fc[faceI]](5,0) -= curSource.z();
                }
            }
        }
        else if
        (
            isA<fixedValueFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pW =
                refCast<fixedValueFvPatchVectorField>
                (
                    W_.boundaryFieldRef()[patchI]
                );

            const vectorField& pWPrev =
                W_.prevIter().boundaryField()[patchI];

            vectorField& pWCorr =
                DW_.boundaryFieldRef()[patchI];

            pWCorr = pW - pWPrev;

            // Diag contribution
            forAll (pW, faceI)
            {
                d[fc[faceI]](3,0) -= pCMQW[faceI].xx()/pDelta[faceI];
                d[fc[faceI]](3,1) -= pCMQW[faceI].xy()/pDelta[faceI];
                d[fc[faceI]](3,2) -= pCMQW[faceI].xz()/pDelta[faceI];

                d[fc[faceI]](4,0) -= pCMQW[faceI].yx()/pDelta[faceI];
                d[fc[faceI]](4,1) -= pCMQW[faceI].yy()/pDelta[faceI];
                d[fc[faceI]](4,2) -= pCMQW[faceI].yz()/pDelta[faceI];

                d[fc[faceI]](5,0) -= pCMQW[faceI].zx()/pDelta[faceI];
                d[fc[faceI]](5,1) -= pCMQW[faceI].zy()/pDelta[faceI];
                d[fc[faceI]](5,2) -= pCMQW[faceI].zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pW, faceI)
            {
                vector curSource =
                (
                    pCMQW[faceI]
                  & (
                        pWCorr[faceI]
                      + vector(SMALL, SMALL, SMALL)
                    )
                )/pDelta[faceI];

                source[fc[faceI]](3,0) -= curSource.x();
                source[fc[faceI]](4,0) -= curSource.y();
                source[fc[faceI]](5,0) -= curSource.z();
            }

            if
            (
                isA<momentBeamRotationNRFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationNRFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationNRFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                scalarField pDelta
                (
                    1.0/mesh().deltaCoeffs().boundaryField()[patchI]
                );
                tensorField invCM(inv(pCMTheta/pDelta + pCMTheta2));

                vectorField A((invCM & (pTheta.moment() - pExplicitM)));
                tensorField B((invCM & (pCMTheta/pDelta)));

                vectorField src((pCMQTheta & A));
                tensorField diag((pCMQTheta & B));

                // Diag contribution
                forAll (pTheta, faceI)
                {
                    d[fc[faceI]](3,3) += diag[faceI].xx();
                    d[fc[faceI]](3,4) += diag[faceI].xy();
                    d[fc[faceI]](3,5) += diag[faceI].xz();

                    d[fc[faceI]](4,3) += diag[faceI].yx();
                    d[fc[faceI]](4,4) += diag[faceI].yy();
                    d[fc[faceI]](4,5) += diag[faceI].yz();

                    d[fc[faceI]](5,3) += diag[faceI].zx();
                    d[fc[faceI]](5,4) += diag[faceI].zy();
                    d[fc[faceI]](5,5) += diag[faceI].zz();
                }

                // Source contribution
                forAll (pTheta, faceI)
                {
                    source[fc[faceI]](3,0) -= src[faceI].x();
                    source[fc[faceI]](4,0) -= src[faceI].y();
                    source[fc[faceI]](5,0) -= src[faceI].z();
                }
            }
            else if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        Theta_.boundaryFieldRef()[patchI]
                    );

                const vectorField& pThetaPrev =
                    Theta_.prevIter().boundaryField()[patchI];

                vectorField& pThetaCorr =
                    DTheta_.boundaryFieldRef()[patchI];

                // pThetaCorrStar = (pTheta - pThetaPrev);
                pThetaCorr = (pLambda & (pTheta - pThetaPrev));

                // Info << patchI << ", " << pThetaCorr << endl;

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector curSource =
                        (pCMQTheta[faceI] & pThetaCorr[faceI]);

                    source[fc[faceI]](3,0) -= curSource.x();
                    source[fc[faceI]](4,0) -= curSource.y();
                    source[fc[faceI]](5,0) -= curSource.z();
                }
            }
        }

        // Add explicit force moment
        forAll (pExplicitMQ, faceI)
        {
            if (true)
            {
                vector correctedOwnExplicitMQ = pExplicitMQ[faceI];

                // correctedOwnExplicitMQ =
                // (
                //     refRMI[fc[faceI]].T()
                //   & (RMI[fc[faceI]].T() & correctedOwnExplicitMQ)
                // );
                // correctedOwnExplicitMQ.x() = 0;
                // correctedOwnExplicitMQ =
                // (
                //     RMI[fc[faceI]]
                //  & (refRMI[fc[faceI]] & correctedOwnExplicitMQ)
                // );

                source[fc[faceI]](3,0) -= correctedOwnExplicitMQ.x();
                source[fc[faceI]](4,0) -= correctedOwnExplicitMQ.y();
                source[fc[faceI]](5,0) -= correctedOwnExplicitMQ.z();
            }
        }
    }
