/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "linearElasticMisesPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticMisesPlastic, 0);
    addToRunTimeSelectionTable
    (
        plasticityStressResultantReturn,
        linearElasticMisesPlastic,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


    
tmp<scalarField>
linearElasticMisesPlastic::yieldFunction
(
    const vectorField& S,
    const scalarField& eqEP
) const
{
    tmp<scalarField> tYieldFunction
    (
        new scalarField(S.size(), 0)
    );

    diagTensor P(1, 3, 3);

    tYieldFunction() = sqrt(S & (P & S))
      - (sigmaY0_.value() + HP_.value()*eqEP);

    return tYieldFunction;
}

scalar linearElasticMisesPlastic::yieldFunction
(
    const vector& S,
    const scalar& eqEP
) const
{
    diagTensor P(1, 3, 3);

    scalar yf = sqrt(S & (P & S))
      - (sigmaY0_.value() + HP_.value()*eqEP); 

    return yf;
}


scalar linearElasticMisesPlastic::yieldFunctionDerivative
(
    const vector& S,
    const scalar& eqEP,
    const scalar& gamma
) const
{
    diagTensor P(1, 3, 3);

    scalar sigmaY = sigmaY0_.value() + HP_.value()*eqEP;

    scalar eqS = vonMisesStress(S);

    vector N = (P & S)/eqS;

    diagTensor bC = barC(eqEP, gamma);

    // scalar EE = beamModel_.E().value();
    // scalar GG = beamModel_.G().value();

    // (
    //     EE/(1.0 + EE*gamma/sigmaY),                    
    //     GG/(1.0 + 3*GG*gamma/sigmaY),
    //     GG/(1.0 + 3*GG*gamma/sigmaY)
    // );

    scalar dFdGamma =
   -(
        (eqS/sigmaY)*(1.0 - gamma*HP_.value()/sigmaY)*(N & (bC & N))
      + HP_.value()
    );

    return dFdGamma;
}


scalar linearElasticMisesPlastic::vonMisesStress
(
    const vector& S
) const
{
    diagTensor P(1, 3, 3);

    scalar eqS = sqrt(S & (P & S));

    return eqS;
}


tmp<scalarField> linearElasticMisesPlastic::vonMisesStress
(
    const vectorField& S
) const
{
    tmp<scalarField> tEqS
    (
        new scalarField(S.size(), 0)
    );
    
    diagTensor P(1, 3, 3);

    tEqS() = sqrt(S & (P & S));

    return tEqS;
}

    
diagTensor linearElasticMisesPlastic::barC
(
    const scalar eqEP,
    const scalar gamma
) const
{
    diagTensor P(1, 3, 3);

    scalar sigmaY = sigmaY0_.value() + HP_.value()*eqEP;

    scalar EE = beamModel_.E().value();
    scalar GG = beamModel_.G().value();

    diagTensor invC(1.0/EE, 1.0/GG, 1.0/GG);

    diagTensor bC = inv(invC + (gamma/sigmaY)*P);

    return bC;
}


tmp<diagTensorField> linearElasticMisesPlastic::barC
(
    const scalarField& eqEP,
    const scalarField& gamma
) const
{
    tmp<diagTensorField> tBarC
    (
        new diagTensorField(eqEP.size(), diagTensor::zero)
    );
    
    diagTensor P(1, 3, 3);

    scalarField sigmaY = sigmaY0_.value() + HP_.value()*eqEP;
                    
    scalar EE = beamModel_.E().value();
    scalar GG = beamModel_.G().value();

    diagTensor invC(1.0/EE, 1.0/GG, 1.0/GG);

    tBarC() = inv(invC + (gamma/sigmaY)*P);

    return tBarC;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
linearElasticMisesPlastic::linearElasticMisesPlastic
(
    const word& name,
    beamModel& beamModel
)
:
    plasticityStressResultantReturn(name, beamModel),
    beamModel_(beamModel),
    eqEP_(beamModel_.mesh().nFaces()),
    EP_(beamModel_.mesh().nFaces()),
    gamma_(beamModel_.mesh().nFaces()),
    plasticN_(beamModel_.mesh().nFaces()),
    sigmaY0_(plasticityProperties().lookup("sigmaY0")),
    HP_(plasticityProperties().lookup("HP")),
    DGammaP_
    (
        IOobject
        (
            "DGammaP",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beamModel.mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    DKP_
    (
        IOobject
        (
            "DKP",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beamModel.mesh(),
        dimensionedVector("0", dimless/dimLength, vector::zero)
    ),
    DQDGamma_
    (
        IOobject
        (
            "DQDGamma",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        beamModel.mesh(),
        dimensionedTensor("0", dimPressure*dimArea, tensor::zero)
    ),
    DMDK_
    (
        IOobject
        (
            "DMDK",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        beamModel.mesh(),
        dimensionedTensor
        (
            "0",
            dimPressure*dimArea*dimArea,
            tensor::zero
        )
    ),
    DMDGamma_
    (
        IOobject
        (
            "DMDGamma",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        beamModel.mesh(),
        dimensionedTensor
        (
            "0",
            dimForce*dimLength,
            tensor::zero
        )
    ),
    DQDK_
    (
        IOobject
        (
            "DQDK",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        beamModel.mesh(),
        dimensionedTensor
        (
            "0",
            dimForce*dimLength,
            tensor::zero
        )
    ),
    nYieldingFaces_(0),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        beamModel.mesh(),
        dimensionedScalar
        (
            "0",
            dimless,
            0
        )
    ),
    writeFaceIDs_(plasticityProperties().lookup("writeFaceIDs")),
    plasticity_
    (
        plasticityProperties().lookupOrDefault<bool>("plasticity", true)
    )
{
    Info << "Creating linearElasticMisesPlastic "
         << "plastic stress resultant return method"
         << endl;

    Info << "HP: " << HP_ << endl;
    
    // Info << "EP:" << EP_ << endl;

    dimensionedTensor CQ
    (
        "CQ",
        beamModel_.EA().dimensions(),
        tensor
        (
            beamModel_.EA().value(), 0, 0,
            0, beamModel_.GA().value(), 0,
            0, 0, beamModel_.GA().value()
        )
    );

    dimensionedTensor CM
    (
        "CM",
        beamModel_.EI().dimensions(),
        tensor
        (
            beamModel_.GJ().value(), 0, 0,
            0, beamModel_.EIyy().value(), 0,
            0, 0, beamModel_.EIzz().value()
        )
    );

    DQDGamma_ = CQ;
    DMDK_ = CM;

    // Read or initialize eq plastic strains
    {
        IFstream eqEPfile
        (
            beamModel_.runTime().timePath()/"eqEP"
        );

        if (eqEPfile.good())
        {
            Info << "Reading eqEP field" << endl;
            eqEPfile >> eqEP_;
        }
        else
        {
            Info << "Initializing eqEP field" << endl;
            forAll(eqEP_, faceI)
            {
                label nLocalPoints =
                    beamModel_.crossSections()[0].points().size();

                eqEP_.set(faceI, scalarField(nLocalPoints, 0));
            }
        }
    }


    // Read or initialize plastic strains
    {
        IFstream EPfile
        (
            beamModel_.runTime().timePath()/"EP"
        );

        if (EPfile.good())
        {
            Info << "Reading EP field" << endl;
            EPfile >> EP_;
        }
        else
        {
            Info << "Initializing EP field" << endl;
            forAll(EP_, faceI)
            {
                label nLocalPoints =
                    beamModel_.crossSections()[0].points().size();

                EP_.set(faceI, vectorField(nLocalPoints, vector::zero));
            }
        }
    }

    
    forAll(eqEP_, faceI)
    {
        label nLocalPoints =
            beamModel_.crossSections()[0].points().size();
        
        // EP_.set(faceI, vectorField(nLocalPoints, vector::zero));
        // eqEP_.set(faceI, scalarField(nLocalPoints, 0));
        
        gamma_.set(faceI, scalarField(nLocalPoints, 0));
        plasticN_.set(faceI, vectorField(nLocalPoints, vector::zero));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

linearElasticMisesPlastic::~linearElasticMisesPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void linearElasticMisesPlastic::correct()
{
    if (timeIndex() != beamModel_.runTime().timeIndex())
    {
        // Info << "Update total plasticity fields" << endl;
        // EP_ += gamma_*plasticN_;
        // eqEP_ += gamma_;

        gamma_ = 0;
        plasticN_ = vector::zero;
        
        activeYield_ =
            dimensionedScalar("0", activeYield_.dimensions(), 0.0);
        
        timeIndex() = beamModel_.runTime().timeIndex();
    }

    const fvMesh& mesh = beamModel_.mesh();
    
    // Look up for fields
    const surfaceVectorField& Gamma =
        beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
        (
            "Gamma"
        );

    const surfaceVectorField& GammaP =
        beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
        (
            "GammaP"
        );

    const surfaceVectorField& K =
        beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
        (
            "K"
        );

    const surfaceVectorField& KP =
        beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
        (
            "KP"
        );

    scalar EE = beamModel_.E().value();
    scalar GG = beamModel_.G().value();
    diagTensor C(EE, GG, GG);
    
    diagTensor P(1, 3, 3);

    dimensionedTensor CQ
    (
        "CQ",
        beamModel_.EA().dimensions(),
        tensor
        (
            beamModel_.EA().value(), 0, 0,
            0, beamModel_.GA().value(), 0,
            0, 0, beamModel_.GA().value()
        )
    );

    dimensionedTensor CM
    (
        "CM",
        beamModel_.EI().dimensions(),
        tensor
        (
            beamModel_.GJ().value(), 0, 0,
            0, beamModel_.EIyy().value(), 0,
            0, 0, beamModel_.EIzz().value()
        )
    );

    const vectorField& GammaPI = GammaP.internalField();
    const vectorField& KPI = KP.internalField();

    const vectorField& oldGammaPI = GammaP.oldTime().internalField();
    const vectorField& oldKPI = KP.oldTime().internalField();

    const unallocLabelList& own = beamModel_.mesh().owner();
    const unallocLabelList& nei = beamModel_.mesh().neighbour();

    const vectorField& GammaI = Gamma.internalField();
    const vectorField& KI = K.internalField();

    nYieldingFaces_ = 0;

    scalar smallF = SMALL;

    forAll(GammaI, faceI)
    {
        vectorField E =
            beamModel_.crossSections()[0]
           .greenLagrangianStrain(GammaI[faceI], KI[faceI]);

        vectorField Etrial = E - EP_[faceI];

        vectorField Strial = (C & Etrial);

        // vector checkQ =
        //     beamModel_.crossSections()[0].resultantForce(Strial);
        // vector checkM =
        //     beamModel_.crossSections()[0].resultantMoment(Strial);

        scalarField Ftrial = yieldFunction(Strial, eqEP_[faceI]);

        bool yielding = false;

        forAll(Ftrial, pointI)
        {
            vector S = Strial[pointI];

            if
            (
                (Ftrial[pointI] > smallF)
             && plasticity_
            )
            {
                // Calculate consistency parameter
                scalar gamma = gamma_[faceI][pointI];
                scalar residual = GREAT;

                label iCorr = 0;

                do
                {
                    scalar eqEP = eqEP_[faceI][pointI] + gamma;

                    diagTensor bC = barC(eqEP, gamma);

                    S = (bC & Etrial[pointI]);

                    scalar F = yieldFunction(S, eqEP);

                    scalar dFdGamma =
                        yieldFunctionDerivative(S, eqEP, gamma);

                    scalar newGamma = gamma - F/dFdGamma;

                    residual = mag(newGamma - gamma)/(mag(newGamma) + SMALL);

                    gamma = newGamma;
                }
                while
                (
                    (residual > 1e-6) && (++iCorr < 100)
                );

                if (iCorr==100)
                {
                    Info << "Num of Newton iterations: " << iCorr
                         << " (" << faceI << ", " << pointI << ")" << endl;
                    Info << residual << ", " << gamma << endl;
                }

                gamma_[faceI][pointI] = gamma;
                plasticN_[faceI][pointI] = (P & S)/vonMisesStress(S);
                yielding = true;
            }
            else
            {
                gamma_[faceI][pointI] = 0;
                plasticN_[faceI][pointI] = (P & S)/vonMisesStress(S);
            }
        }

        if (yielding)
        {
            // Info << QI[faceI] << ", " << checkQ << endl;
            // Info << MI[faceI] << ", " << checkM << endl;

            // Calculate resultant plastic strain correction
            scalarField eqEpl = eqEP_[faceI] + gamma_[faceI];
            vectorField Epl = EP_[faceI] + gamma_[faceI]*plasticN_[faceI];
            vectorField Eel = E - Epl;
            vectorField S = (C & Eel);

            vector newQ = beamModel_.crossSections()[0].resultantForce(S);
            vector newM = beamModel_.crossSections()[0].resultantMoment(S);

            vector newGammaP = GammaI[faceI] - (inv(CQ.value()) & newQ);
            vector newKP = KI[faceI] - (inv(CM.value()) & newM);

            DGammaP_[faceI] = newGammaP - GammaPI[faceI];
            DKP_[faceI] = newKP - KPI[faceI];

            // Calculate plastic modulus
            diagTensorField bC = barC(eqEpl, gamma_[faceI]);
            const vectorField& N = plasticN_[faceI];

            scalarField sigmaY =
                sigmaY0_.value() + HP_.value()*eqEpl;

            scalarField beta =
                HP_.value()
               /(1.0 - (gamma_[faceI]/sigmaY)*HP_.value());

            diagTensorField dSdE
            (
                Ftrial.size(),
                diagTensor(EE, GG, GG)
            );

            // forAll(Ftrial, pI)
            // {
            //     if (Ftrial[pI] > smallF)
            //     {
            //         dSdE[pI] = bC[pI]
            //           - (((bC[pI] & N[pI]) & N[pI]) * bC[pI])
            //            /((N[pI] & (bC[pI] & N[pI])) + beta[pI] + SMALL);
            //     }
            // }

            beamModel_.crossSections()[0].resultantTangentMatrices
            (
                dSdE,
                DQDGamma_[faceI],
                DQDK_[faceI],
                DMDGamma_[faceI],
                DMDK_[faceI]
            );

            // Info << faceI << endl;
            // Info << DQDGamma_[faceI] << ", " <<  CQ.value() << endl;
            // Info << DMDK_[faceI] << ", " <<  CM.value() << endl;

            activeYield_[own[faceI]] = 1;
            activeYield_[nei[faceI]] = 1;

            nYieldingFaces_++;

            // DGammaP_[faceI] = vector::zero;
            // DKP_[faceI] = vector::zero;

            // DQDGamma_[faceI] = CQ.value();
            // DMDK_[faceI] = CM.value();

            // DMDGamma_[faceI] = tensor::zero;
            // DQDK_[faceI] = tensor::zero;
        }
        else
        {
            // // Calculate resultant plastic strain correction
            // scalarField eqEpl = eqEP_[faceI];
            // vectorField Epl = EP_[faceI];
            // vectorField Eel = E - Epl;
            // vectorField S = (C & Eel);

            // vector newQ = beamModel_.crossSections()[0].resultantForce(S);
            // vector newM = beamModel_.crossSections()[0].resultantMoment(S);

            // vector newGammaP = GammaI[faceI] - (inv(CQ.value()) & newQ);
            // vector newKP = KI[faceI] - (inv(CM.value()) & newM);

            // DGammaP_[faceI] = newGammaP - GammaPI[faceI];
            // DKP_[faceI] = newKP - KPI[faceI];

            DGammaP_[faceI] = oldGammaPI[faceI] - GammaPI[faceI];
            DKP_[faceI] = oldKPI[faceI] - KPI[faceI];

            DQDGamma_[faceI] = CQ.value();
            DMDK_[faceI] = CM.value();

            DMDGamma_[faceI] = tensor::zero;
            DQDK_[faceI] = tensor::zero;
        }
    }

    forAll(Gamma.boundaryField(), patchI)
    {
        const vectorField& pGammaI = Gamma.boundaryField()[patchI];
        const vectorField& pKI = K.boundaryField()[patchI];

        const vectorField& pGammaPI = GammaP.boundaryField()[patchI];
        const vectorField& pKPI = KP.boundaryField()[patchI];

        const vectorField& pOldGammaPI =
            GammaP.oldTime().boundaryField()[patchI];
        const vectorField& pOldKPI =
            KP.oldTime().boundaryField()[patchI];

        const labelList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(pGammaI, pFaceI)
        {
            label faceI =
                beamModel_.mesh().boundaryMesh()[patchI].start()
              + pFaceI;

            vectorField E =
                beamModel_.crossSections()[0]
               .greenLagrangianStrain(pGammaI[pFaceI], pKI[pFaceI]);

            vectorField Etrial = E - EP_[faceI];

            vectorField Strial = (C & Etrial);

            scalarField Ftrial = yieldFunction(Strial, eqEP_[faceI]);

            bool yielding = false;

            forAll(Ftrial, pointI)
            {
                vector S = Strial[pointI];

                if (Ftrial[pointI] > smallF)
                {
                    // Calculate consistency parameter
                    scalar gamma = gamma_[faceI][pointI];
                    scalar residual = GREAT;

                    label iCorr = 0;

                    do
                    {
                        scalar eqEP = eqEP_[faceI][pointI] + gamma;

                        diagTensor bC = barC(eqEP, gamma);

                        S = (bC & Etrial[pointI]);

                        scalar F = yieldFunction(S, eqEP);

                        scalar dFdGamma =
                            yieldFunctionDerivative(S, eqEP, gamma);

                        scalar newGamma = gamma - F/dFdGamma;

                        residual =
                            mag(newGamma - gamma)/(mag(newGamma) + SMALL);

                        gamma = newGamma;
                    }
                    while
                    (
                        (residual > 1e-6) && (++iCorr < 100)
                    );
                
                    gamma_[faceI][pointI] = gamma;
                    plasticN_[faceI][pointI] = (P & S)/vonMisesStress(S);
                    yielding = true;
                }
                else
                {
                    gamma_[faceI][pointI] = 0;
                    plasticN_[faceI][pointI] = (P & S)/vonMisesStress(S);
                }
            }

            if (yielding)
            {
                // Calculate resultant plastic strain correction
                scalarField eqEpl = eqEP_[faceI] + gamma_[faceI];
                vectorField Epl = EP_[faceI] + gamma_[faceI]*plasticN_[faceI];
                vectorField Eel = E - Epl;
                vectorField S = (C & Eel);

                vector newQ = beamModel_.crossSections()[0].resultantForce(S);
                vector newM = beamModel_.crossSections()[0].resultantMoment(S);

                vector newGammaP = pGammaI[pFaceI] - (inv(CQ.value()) & newQ);
                vector newKP = pKI[pFaceI] - (inv(CM.value()) & newM);

                DGammaP_.boundaryField()[patchI][pFaceI] =
                    newGammaP - pGammaPI[pFaceI];
                DKP_.boundaryField()[patchI][pFaceI] =
                    newKP - pKPI[pFaceI];

                // // Calculate plastic modulus
                // diagTensorField bC = barC(eqEpl, gamma_[faceI]);
                // const vectorField& N = plasticN_[faceI];

                // scalarField sigmaY =
                //     sigmaY0_.value() + HP_.value()*eqEpl;
                
                // scalarField beta =
                //     HP_.value()
                //    /(1.0 + (gamma_[faceI]/sigmaY)*HP_.value());
                
                // tensorField dSdE =
                //     bC - (bC & N)*(N & bC)/((N & (bC & N)) + beta + SMALL);

                // beamModel_.crossSections()[0].resultantTangentMatrices
                // (
                //     dSdE,
                //     DQDGamma_[faceI],
                //     DQDK_[faceI],
                //     DMDGamma_[faceI],
                //     DMDK_[faceI]
                // );

                // DGammaP_[faceI] = vector::zero;
                // DKP_[faceI] = vector::zero;

                DQDGamma_.boundaryField()[patchI][pFaceI] = CQ.value();
                DMDK_.boundaryField()[patchI][pFaceI] = CM.value();

                DMDGamma_.boundaryField()[patchI][pFaceI] = tensor::zero;
                DQDK_.boundaryField()[patchI][pFaceI] = tensor::zero;

                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const processorFvPatch& procPatch =
                        refCast<const processorFvPatch>
                        (
                            mesh.boundary()[patchI]
                        );
                    if (procPatch.master())
                    {
                        nYieldingFaces_++;
                    }
                }
                else
                {
                    activeYield_.boundaryField()[patchI][pFaceI] = 1;
                    nYieldingFaces_++;
                }
                activeYield_[faceCells[pFaceI]] = 1;
            }
            else
            {
                // // Calculate resultant plastic strain correction
                // scalarField eqEpl = eqEP_[faceI];
                // vectorField Epl = EP_[faceI];
                // vectorField Eel = E - Epl;
                // vectorField S = (C & Eel);

                // vector newQ = beamModel_.crossSections()[0].resultantForce(S);
                // vector newM = beamModel_.crossSections()[0].resultantMoment(S);

                // vector newGammaP = pGammaI[pFaceI] - (inv(CQ.value()) & newQ);
                // vector newKP = pKI[pFaceI] - (inv(CM.value()) & newM);

                // DGammaP_.boundaryField()[patchI][pFaceI] =
                //     newGammaP - pGammaPI[pFaceI];
                // DKP_.boundaryField()[patchI][pFaceI] =
                //     newKP - pKPI[pFaceI];
                
                DGammaP_.boundaryField()[patchI][pFaceI] =
                    pOldGammaPI[pFaceI] - pGammaPI[pFaceI];
                DKP_.boundaryField()[patchI][pFaceI] =
                    pOldKPI[pFaceI] - pKPI[pFaceI];

                DQDGamma_.boundaryField()[patchI][pFaceI] = CQ.value();
                DMDK_.boundaryField()[patchI][pFaceI] = CM.value();

                DMDGamma_.boundaryField()[patchI][pFaceI] = tensor::zero;
                DQDK_.boundaryField()[patchI][pFaceI] = tensor::zero;
            }
        }
    }

    reduce(nYieldingFaces_, sumOp<label>());

    // Info << "nPlasticFaces: " << nYieldingFaces_ << endl;
}

void linearElasticMisesPlastic::updateYieldStress()
{
    Info << "Update total plasticity fields" << endl;

    EP_ += gamma_*plasticN_;
    eqEP_ += gamma_;

    // Reset plasticity fields
    DGammaP_ =
        dimensionedVector("zero", DGammaP_.dimensions(), vector::zero);
    DKP_ =
        dimensionedVector("zero", DKP_.dimensions(), vector::zero);

    dimensionedTensor CQ
    (
        "CQ",
        beamModel_.EA().dimensions(),
        tensor
        (
            beamModel_.EA().value(), 0, 0,
            0, beamModel_.GA().value(), 0,
            0, 0, beamModel_.GA().value()
        )
    );

    dimensionedTensor CM
    (
        "CM",
        beamModel_.EI().dimensions(),
        tensor
        (
            beamModel_.GJ().value(), 0, 0,
            0, beamModel_.EIyy().value(), 0,
            0, 0, beamModel_.EIzz().value()
        )
    );

    DQDGamma_ = CQ;
    DMDK_ = CM;

    DQDK_ = dimensionedTensor("zero", dimForce*dimLength, tensor::zero);
    DMDGamma_ = dimensionedTensor("zero", dimForce*dimLength, tensor::zero);
}

void linearElasticMisesPlastic::writeInfo() const
{
    Info << "Number of yielding faces: "
         << nYieldingFaces_ << endl;
}

void linearElasticMisesPlastic::writeFields() const
{
    Info << "Write total plasticity fields" << endl;
    
    // dimensionedTensor CQ
    // (
    //     "CQ",
    //     beamModel_.EA().dimensions(),
    //     tensor
    //     (
    //         beamModel_.EA().value(), 0, 0,
    //         0, beamModel_.GA().value(), 0,
    //         0, 0, beamModel_.GA().value()
    //     )
    // );

    // dimensionedTensor CM
    // (
    //     "CM",
    //     beamModel_.EI().dimensions(),
    //     tensor
    //     (
    //         beamModel_.GJ().value(), 0, 0,
    //         0, beamModel_.EIyy().value(), 0,
    //         0, 0, beamModel_.EIzz().value()
    //     )
    // );
    
    const surfaceVectorField& Gamma =
        beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
        (
            "Gamma"
        );

    const surfaceVectorField& K =
        beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
        (
            "K"
        );

    const vectorField& GammaI = Gamma.internalField();
    const vectorField& KI = K.internalField();
    
    // const surfaceVectorField& Q =
    //     beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
    //     (
    //         "Q"
    //     );

    // const surfaceVectorField& M =
    //     beamModel_.mesh().objectRegistry::lookupObject<surfaceVectorField>
    //     (
    //         "M"
    //     );

    // const vectorField& QI = Q.internalField();
    // const vectorField& MI = M.internalField();

    // Write cross-section fields
    forAll(writeFaceIDs_, faceI)
    {
        label curFace = writeFaceIDs_[faceI];

        // Create directory if does not exist.
        fileName vtkDir(beamModel_.runTime().path()/"VTK");
        mkDir(vtkDir);

        fileName plasticityDir
        (
            vtkDir/"plasticity"
        );
        mkDir(plasticityDir);

        OStringStream FileName;
        FileName() << "face-" << curFace << "_"
                   << beamModel_.runTime().timeIndex() << ".vtk";
        fileName vtkFileName(plasticityDir/word(FileName.str()));

        // Calculate resultant plastic strain correction
        scalarField eqEpl = eqEP_[curFace]; // It is already updated
        // scalarField eqEpl = eqEP_[curFace] + gamma_[curFace];
        vectorField Epl = EP_[curFace]; // It is already updated
        // vectorField Epl = EP_[curFace] + gamma_[curFace]*plasticN_[curFace];

        vectorField E =
            beamModel_.crossSections()[0]
           .greenLagrangianStrain(GammaI[curFace], KI[curFace]);

        vectorField Eel = E - Epl;

        scalar EE = beamModel_.E().value();
        scalar GG = beamModel_.G().value();
        diagTensor C(EE, GG, GG);

        vectorField S = (C & Eel);

        // vector Q = (CQ.value() & GammaI[curFace]);
        // vector M = (CM.value() & KI[curFace]);
        
        // vector intQ = beamModel_.crossSections()[0].resultantForce(S);
        // vector intM = beamModel_.crossSections()[0].resultantMoment(S);

        // Info << Q << ", " << intQ << endl;
        // Info << M << ", " << intM << endl;
       
        scalarField eqS = vonMisesStress(S);

        beamModel_.crossSections()[0].writeVTK
        (
            vtkFileName,
            eqEpl,
            eqS,
            S
        );

        OStringStream FileNameVtp;
        FileNameVtp() << "face-" << curFace
                      << "_" << beamModel_.runTime().timeIndex()
                      << ".vtp";

        fileName vtpFileName
        (
            word(FileNameVtp.str())
        );

        // Write temporal collection data
        {
            OStringStream faceNpvd;
            faceNpvd() << "face-" << curFace << ".pvd";

            fileName collectionFileName
            (
                plasticityDir/faceNpvd.str()
            );
            ifstream collectionFile(collectionFileName);
            // IFstream collectionFile(collectionFileName);

            fileName newCollectionFileName
            (
                plasticityDir/"new.pvd"
            );

            if (collectionFile.good())
            {
                ofstream newCollectionFile
                (
                    newCollectionFileName
                );

                // Add to existing collectin file
                label lineIndex = 0;
                do
                {
                    lineIndex++;

                    std::string line;
                    std::getline(collectionFile, line);

                    if
                    (
                        (lineIndex < 3)
                     || (
                            line.find("DataSet")
                         != std::string::npos
                        )
                    )
                    {
                        newCollectionFile << line << '\n';
                    }
                }
                while(!collectionFile.eof());

                // Add current pulley data
                newCollectionFile
                    << "    <DataSet timestep=\""
                    << beamModel_.runTime().value()
                    << "\" group=\"\" part=\"0\" file=\""
                    << vtpFileName << "\"/>" << '\n';

                // Add last two lines
                newCollectionFile << "  </Collection>" << '\n';
                newCollectionFile << "</VTKFile>";
            }
            else
            {
                OFstream newCollectionFile
                (
                    newCollectionFileName
                );
                        
                // Add first two lines
                newCollectionFile
                    << "<VTKFile type=\"Collection\" version=\"0.1\" "
                    << "byte_order=\"LittleEndian\">" << endl;
                newCollectionFile << "  <Collection>" << endl;

                // Add current pulley data
                newCollectionFile
                    << "    <DataSet timestep=\""
                    << beamModel_.runTime().value()
                    << "\" group=\"\" part=\"0\" file="
                    << vtpFileName << "/>" << endl;
                
                // Add last two lines
                newCollectionFile << "  </Collection>" << endl;
                newCollectionFile << "</VTKFile>";                            
            }

            mv(newCollectionFileName, collectionFileName);
        }

        // Write traction profile
        {
            OStringStream faceData;
            faceData() << "crossSectionTraction_" << curFace << ".dat";

            fileName fName
            (
                beamModel_.runTime().timePath()/faceData.str()
            );
            OFstream file(fName);

            vectorField csPoints =
                beamModel_.crossSections()[0].points();

            forAll(csPoints, pointI)
            {
                file << csPoints[pointI].y() << " "
                     << S[pointI].x() << " "
                     << S[pointI].y() << " "
                     << S[pointI].z() << " "
                     << eqEpl[pointI] << endl;
            }

            // file.closed();
        }
    }

    // Write cross-section eq plastic streains
    {
        OFstream eqEPfile
        (
            beamModel_.runTime().timePath()/"eqEP"
        );

        eqEPfile << eqEP_ << endl;
    }

    // Write cross-section plastic streains
    {
        OFstream EPfile
        (
            beamModel_.runTime().timePath()/"EP"
        );

        EPfile << EP_ << endl;
    }
}

} // end of namespace

// ************************************************************************* //
