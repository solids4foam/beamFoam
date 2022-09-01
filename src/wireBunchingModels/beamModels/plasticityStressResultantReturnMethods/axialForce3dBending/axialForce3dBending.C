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

#include "axialForce3dBending.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(axialForce3dBending, 0);
    addToRunTimeSelectionTable
    (
        plasticityStressResultantReturn,
        axialForce3dBending,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar axialForce3dBending::yieldFunction
(
    const scalar N,
    const vector& M,
    const scalar NP
)
{
    return sqr(N/NP) + CNM_.value()*sqrt(sqr(M.y())+sqr(M.z()))/NP - 1.0;
}

void axialForce3dBending::calcPlasticStrainIncrement
(
    scalar& DCP,
    scalar& DEP,
    vector& DEpsilonP,
    tensor& DSDEpsilon,
    const scalar N,
    const vector& M,
    const scalar NP,
    const scalar CP,
    const vector dEpsilonP,
    const scalar dEP,
    const scalar Ftrial
)
{
    vector DFDS
    (
        2*N/sqr(NP),
        CNM_.value()*(M.y()/sqrt(sqr(M.y())+sqr(M.z())))/NP,
        CNM_.value()*(M.z()/sqrt(sqr(M.y())+sqr(M.z())))/NP
    );

    tensor D2FDS2
    (
        2.0/sqr(NP), 0, 0,
        0, CNM_.value()*sqr(M.z()/(sqr(M.y())+sqr(M.z())))*sqrt(sqr(M.y())+sqr(M.z()))/NP,
       -CNM_.value()*(M.y()*M.z()/sqr(sqr(M.y())+sqr(M.z())))*sqrt(sqr(M.y())+sqr(M.z()))/NP,
        0, -CNM_.value()*(M.y()*M.z()/sqr(sqr(M.y())+sqr(M.z())))*sqrt(sqr(M.y())+sqr(M.z()))/NP,
        CNM_.value()*sqr(M.y()/(sqr(M.y())+sqr(M.z())))*sqrt(sqr(M.y())+sqr(M.z()))/NP
    );

    vector D2FDSDNP
    (
        -4*N/pow(NP,3),
        -CNM_.value()*(M.y()/sqrt(sqr(M.y())+sqr(M.z())))/sqr(NP),
        -CNM_.value()*(M.z()/sqrt(sqr(M.y())+sqr(M.z())))/sqr(NP)
    );

    scalar DFDNP =
  - (
        CNM_.value()*sqrt(sqr(M.y())+sqr(M.z()))/sqr(NP)
      + 2*sqr(N)/pow(NP,3)
    );

    scalar D2FDNP2 =
    (
        2*CNM_.value()*sqrt(sqr(M.y())+sqr(M.z()))/pow(NP,3)
      + 6*sqr(N)/pow(NP,4)
    );

    vector RS = CP*DFDS - dEpsilonP;
    scalar RNP = -CP*DFDNP - dEP;

    label m = 1;
    if (dEP < 0)
    {
        m = -1;
    }

    const beamModel& bm = beamModel_;

    tensor invD
    (
        1.0/bm.EA().value(), 0, 0,
        0, 1.0/bm.EIyy().value(), 0,
        0, 0, 1.0/bm.EIzz().value()
    );

    tensor T0 = inv(invD + CP*D2FDS2);

    scalar alpha = 1.0
      + (EP_.value()*bm.A().value()/m)
       *(
            CP*D2FDNP2
          - sqr(CP)*(D2FDSDNP & (T0 & D2FDSDNP))
        );

    scalar M11 = EP_.value()*bm.A().value()/(m*alpha);

    vector M01 = - CP * (D2FDSDNP & T0) * M11;
    
    tensor M00 = (T0 & (tensor::I - CP*(D2FDSDNP*M01)));

    scalar B =
        (DFDS & ((M00 & RS) - (M01*RNP)))
      + (DFDNP*((M01 & RS) - (M11*RNP)));

    scalar C =
        (DFDS & ((M00 & DFDS) + (M01*DFDNP)))
      + (DFDNP*((M01 & DFDS) + (M11*DFDNP)));

    DCP = (Ftrial - B)/C;

    vector DS =
    (
        (M00 & (-RS - DCP*DFDS))
      + (M01 * (RNP - DCP*DFDNP))
    );

    scalar DNP =
    (
        (M01 & (-RS - DCP*DFDS))
      + (M11 * (RNP - DCP*DFDNP))
    );

    DEpsilonP = -(invD & DS);
    DEP = (m/(bm.A().value()*EP_.value()+SMALL))*DNP;

    // Calculate elasto-plastic moduli
    {
        vector D = (DFDS & M00) + DFDNP*M01;
        
        scalar E =
            (DFDS & (M00 & DFDS))
          + 2*(DFDS & M01)*DFDNP
          + sqr(DFDNP)*M11;
        
        vector beta = D/E;

        DSDEpsilon =
            (M00 & (tensor::I - (DFDS*beta)))
          - DFDNP*(M01*beta);
    }
}


void axialForce3dBending::calcPerfectPlasticStrainIncrement
(
    scalar& DCP,
    scalar& DEP,
    vector& DEpsilonP,
    tensor& DSDEpsilon,
    const scalar N,
    const vector& M,
    const scalar NP,
    const scalar CP,
    const vector dEpsilonP,
    const scalar dEP,
    const scalar Ftrial
)
{
    const beamModel& bm = beamModel_;
        
    tensor invD
    (
        1.0/bm.EA().value(), 0, 0,
        0, 1.0/bm.EIyy().value(), 0,
        0, 0, 1.0/bm.EIzz().value()
    );

    vector DS = vector::zero;

    scalar newN = N;
    vector newM = M;
    scalar newFtrial = Ftrial;
    scalar newCP = CP;
    vector newdEpsilonP = dEpsilonP;

    label iCorr = 0;
    do
    {
        vector DFDS
        (
            2*newN/sqr(NP),
            CNM_.value()*(newM.y()/sqrt(sqr(newM.y())+sqr(newM.z())))/NP,
            CNM_.value()*(newM.z()/sqrt(sqr(newM.y())+sqr(newM.z())))/NP
        );

        tensor D2FDS2
        (
            2.0/sqr(NP), 0, 0,
            0, CNM_.value()*sqr(newM.z()/(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP,
           -CNM_.value()*(newM.y()*M.z()/sqr(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP,
            0, -CNM_.value()*(newM.y()*newM.z()/sqr(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP,
            CNM_.value()*sqr(newM.y()/(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP
        );

        vector RS = newCP*DFDS - newdEpsilonP;

        tensor M0 = inv(invD + newCP*D2FDS2);

        scalar B = (DFDS & (M0 & RS));

        scalar C = (DFDS & (M0 & DFDS));

        DCP = (newFtrial - B)/C;

        DS = -(M0 & (RS + DCP*DFDS));

        DEpsilonP = -(invD & DS);
        DEP = 0.0;

        newCP += DCP;
        newN += DS.x();
        newM.y() += DS.y();
        newM.z() += DS.z();
        newFtrial = yieldFunction(newN, newM, NP);
        newdEpsilonP += DEpsilonP;
    }
    while
    (
        false
     //    mag(newFtrial)>SMALL
     // && ++iCorr < 100
    );

    // if (iCorr >= 100)
    // {
    //     Info << newFtrial << endl;
    //     Info << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    // }

    DCP = newCP - CP;
    DEP = 0;
    DEpsilonP = newdEpsilonP - dEpsilonP;
        
    // Calculate elasto-plastic moduli
    {
        scalar newN = N;
        vector newM = M;
        scalar newCP = CP;
    
        // vector newS(newN, newM.y(), newM.z());        
        // vector newM(M.x(), newS.y(), newS.z());
        // scalar newN = newS.x();
        // scalar newCP = CP + DCP;

        
        
        vector newDFDS
        (
            2*newN/sqr(NP),
            CNM_.value()*(newM.y()/sqrt(sqr(newM.y())+sqr(newM.z())))/NP,
            CNM_.value()*(newM.z()/sqrt(sqr(newM.y())+sqr(newM.z())))/NP
        );
    
        tensor newD2FDS2
        (
            2.0/sqr(NP), 0, 0,
            0, CNM_.value()*sqr(newM.z()/(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP,
            -CNM_.value()*(newM.y()*newM.z()/sqr(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP,
            0, -CNM_.value()*(newM.y()*newM.z()/sqr(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP,
            CNM_.value()*sqr(newM.y()/(sqr(newM.y())+sqr(newM.z())))*sqrt(sqr(newM.y())+sqr(newM.z()))/NP
        );
    
        tensor newM0 = inv(invD + newCP*newD2FDS2);

        vector D = (newDFDS & newM0);

        scalar E = (newDFDS & (newM0 & newDFDS));

        vector beta = D/E;

        DSDEpsilon = (newM0 & (tensor::I - (newDFDS*beta)));
    }
}

    
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
axialForce3dBending::axialForce3dBending
(
    const word& name,
    beamModel& beamModel
)
:
    plasticityStressResultantReturn(name, beamModel),
    beamModel_(beamModel),
    eP_
    (
        IOobject
        (
            "eP",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        beamModel.mesh(),
        dimensionedScalar("0", dimless, 0)
    ),
    CP_
    (
        IOobject
        (
            "CP",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        beamModel.mesh(),
        dimensionedScalar("0", dimless, 0)
    ),
    sigmaY0_(plasticityProperties().lookup("sigmaY0")),
    EP_(plasticityProperties().lookup("Ep")),
    CNM_("CNM", dimless/dimLength, (3*M_PI) / (4*beamModel.R()) ),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            beamModel.runTime().timeName(),
            beamModel.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beamModel.mesh(),
        sigmaY0_
    ),
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
    )
{
    Info << "Creating axialForce3dBending "
         << "plastic stress resultant return method"
         << endl;

    // Info << "EP:" << EP_ << endl;

    if (beamModel_.crossSections().size())
    {
        CNM_.value() = beamModel_.crossSections()[0].CNM();
    }
    
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

axialForce3dBending::~axialForce3dBending()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void axialForce3dBending::correct()
{
    if (timeIndex() != beamModel_.runTime().timeIndex())
    {
        CP_ = dimensionedScalar("0", CP_.dimensions(), 0.0);

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

    scalarField& CPI = CP_.internalField();

    scalarField& ePI = eP_.internalField();
    const scalarField& oldEPI = eP_.oldTime().internalField();

    surfaceVectorField Qtrial = (CQ & (Gamma - GammaP.oldTime()));
    const vectorField& QtrialI = Qtrial.internalField();

    surfaceVectorField Q = (CQ & (Gamma - GammaP));
    const vectorField& QI = Q.internalField();

    surfaceVectorField Mtrial = (CM & (K - KP.oldTime()));
    const vectorField& MtrialI = Mtrial.internalField();

    surfaceVectorField M = (CM & (K - KP));
    const vectorField& MI = M.internalField();
    
    const vectorField& GammaPI = GammaP.internalField();
    const vectorField& KPI = KP.internalField();

    const vectorField& oldGammaPI = GammaP.oldTime().internalField();
    const vectorField& oldKPI = KP.oldTime().internalField();

    const unallocLabelList& own = beamModel_.mesh().owner();
    const unallocLabelList& nei = beamModel_.mesh().neighbour();
    
    nYieldingFaces_ = 0;

    scalar smallF = SMALL;

    bool pureBending
    (
        plasticityProperties().lookupOrDefault<bool>("pureBending", false)
    );

    forAll(CPI, faceI)
    {
        scalar Ntrial = QtrialI[faceI].x();
        if (pureBending)
        {
            Ntrial = 0;
        }
        
        const vector& Mtrial = MtrialI[faceI];

        // Calculate current full plastic force
        scalar NP =
            beamModel_.A().value()
           *(sigmaY0_.value() + EP_.value()*mag(ePI[faceI]));

        scalar Ftrial = yieldFunction(Ntrial, Mtrial, NP);

        if (Ftrial > smallF)
        {
            scalar N = QI[faceI].x();
            if (pureBending)
            {
                N = 0;
            }
            const vector& M = MI[faceI];

            scalar F = yieldFunction(N, M, NP);

            vector dEpsilonP
            (
                GammaPI[faceI].x() - oldGammaPI[faceI].x(),
                KPI[faceI].y() - oldKPI[faceI].y(),
                KPI[faceI].z() - oldKPI[faceI].z()
            );
            scalar dEP(ePI[faceI] - oldEPI[faceI]);

            // Correctors
            scalar DCP = 0;
            scalar DEP = 0;
            vector DEpsilonP = vector::zero;
            tensor DSDEpsilon = tensor::zero;

            const scalar& CP = CPI[faceI];

            if (mag(EP_.value()) > SMALL)
            {
                calcPlasticStrainIncrement
                (
                    DCP,
                    DEP,
                    DEpsilonP,
                    DSDEpsilon,
                    N,
                    M,
                    NP,
                    CP,
                    dEpsilonP,
                    dEP,
                    F
                );
            }
            else
            {
                // Info << "test----------------------------" << endl;
                calcPerfectPlasticStrainIncrement
                (
                    DCP,
                    DEP,
                    DEpsilonP,
                    DSDEpsilon,
                    N,
                    M,
                    NP,
                    CP,
                    dEpsilonP,
                    dEP,
                    F
                );
            }

            CPI[faceI] += DCP;
            ePI[faceI] += DEP;

            DGammaP_[faceI].x() = DEpsilonP.x();
            DGammaP_[faceI].y() = 0;
            DGammaP_[faceI].z() = 0;
            DKP_[faceI].x() = 0;
            DKP_[faceI].y() = DEpsilonP.y();
            DKP_[faceI].z() = DEpsilonP.z();

            DQDGamma_[faceI] = //CQ.value();
                tensor
                (
                    DSDEpsilon.xx(), 0, 0,
                    0, CQ.value().yy(), 0,
                    0, 0, CQ.value().zz()
                );

            DMDK_[faceI] = //CM.value();
                tensor
                (
                    CM.value().xx(), 0,               0,
                    0,               DSDEpsilon.yy(), DSDEpsilon.yz(),
                    0,               DSDEpsilon.zy(), DSDEpsilon.zz()
                );

            DMDGamma_[faceI] = //tensor::zero;
                tensor
                (
                    0,               0, 0,
                    DSDEpsilon.yx(), 0, 0,
                    DSDEpsilon.zx(), 0, 0
                );

            DQDK_[faceI] = //tensor::zero;
                tensor
                (
                    0, DSDEpsilon.xy(), DSDEpsilon.xz(),
                    0, 0, 0,
                    0, 0, 0
                );

            // if (faceI == 0)
            // {
            //     Info << "Ftrial: " << Ftrial << endl;
            //     Info << "F: " << F << endl;
            //     Info << "NP: " << NP << endl;
            //     Info << "N: " << N << endl;
            //     Info << "M: " << M << endl;
            //     Info << "DCP: " << DCP << endl;
            //     Info << "DEP: " << DEP << endl;
            //     Info << "dEP: " << dEP << endl;
            //     Info << "DEpsilonP: " << DEpsilonP << endl;
            //     Info << "dEpsilonP: " << dEpsilonP << endl;
            //     Info << "DSDEpsilon: " << DSDEpsilon << endl;
            //     Info << "CQ: " << CQ.value() << endl;
            //     Info << "CM: " << CM.value() << endl;
            //     Info << "DQDGamma_[faceI]: " << DQDGamma_[faceI] << endl;
            //     Info << "DMDK_[faceI]: " << DMDK_[faceI] << endl;
            // }

            nYieldingFaces_++;
            activeYield_[own[faceI]] = 1;
            activeYield_[nei[faceI]] = 1;
        }
        else
        {
            // DGammaP_[faceI] = vector::zero;
            // DKP_[faceI] = vector::zero;
            DGammaP_[faceI] = oldGammaPI[faceI] - GammaPI[faceI];
            DKP_[faceI] = (oldKPI[faceI] - KPI[faceI]);

            DQDGamma_[faceI] = CQ.value();
            DMDK_[faceI] = CM.value();
            
            DMDGamma_[faceI] = tensor::zero;
            DQDK_[faceI] = tensor::zero;
        }
    }

    forAll(CP_.boundaryField(), patchI)
    {
        scalarField& CPI = CP_.boundaryField()[patchI];

        scalarField& ePI = eP_.boundaryField()[patchI];
        const scalarField& oldEPI = eP_.oldTime().boundaryField()[patchI];

        const vectorField& QtrialI = Qtrial.boundaryField()[patchI];
        const vectorField& MtrialI = Mtrial.boundaryField()[patchI];
        
        const vectorField& QI = Q.boundaryField()[patchI];
        const vectorField& MI = M.boundaryField()[patchI];

        const vectorField& GammaPI = GammaP.boundaryField()[patchI];
        const vectorField& KPI = KP.boundaryField()[patchI];

        const vectorField& oldGammaPI =
            GammaP.oldTime().boundaryField()[patchI];
        const vectorField& oldKPI =
            KP.oldTime().boundaryField()[patchI];

        const labelList& faceCells =
            beamModel_.mesh().boundary()[patchI].faceCells();
        
        forAll(CPI, faceI)
        {
            scalar Ntrial = QtrialI[faceI].x();
            if (pureBending)
            {
                Ntrial = 0;
            }

            const vector& Mtrial = MtrialI[faceI];

            // Calculate current full plastic force
            scalar NP =
                beamModel_.A().value()
               *(sigmaY0_.value() + EP_.value()*mag(ePI[faceI]));

            scalar Ftrial = yieldFunction(Ntrial, Mtrial, NP);

            if (Ftrial > smallF)
            {
                scalar N = QI[faceI].x();
                if (pureBending)
                {
                    N = 0;
                }
            
                const vector& M = MI[faceI];

                scalar F = yieldFunction(N, M, NP);
            
                vector dEpsilonP
                (
                    GammaPI[faceI].x() - oldGammaPI[faceI].x(),
                    KPI[faceI].y() - oldKPI[faceI].y(),
                    KPI[faceI].z() - oldKPI[faceI].z()
                );
                scalar dEP(ePI[faceI] - oldEPI[faceI]);

                scalar DCP = 0;
                scalar DEP = 0;
                vector DEpsilonP = vector::zero;
                tensor DSDEpsilon = tensor::zero;
                // diagTensor DSDEpsilon = diagTensor::zero;

                const scalar& CP = CPI[faceI];

                if (mag(EP_.value()) > SMALL)
                {
                    calcPlasticStrainIncrement
                    (
                        DCP,
                        DEP,
                        DEpsilonP,
                        DSDEpsilon,
                        N,
                        M,
                        NP,
                        CP,
                        dEpsilonP,
                        dEP,
                        F
                    );
                }
                else
                {
                    // Info << "btest------------------------" << endl;
                    calcPerfectPlasticStrainIncrement
                    (
                        DCP,
                        DEP,
                        DEpsilonP,
                        DSDEpsilon,
                        N,
                        M,
                        NP,
                        CP,
                        dEpsilonP,
                        dEP,
                        F
                    );
                }

                CPI[faceI] += DCP;
                ePI[faceI] += DEP;

                DGammaP_.boundaryField()[patchI][faceI].x() = DEpsilonP.x();

                DKP_.boundaryField()[patchI][faceI].y() = DEpsilonP.y();
                DKP_.boundaryField()[patchI][faceI].z() = DEpsilonP.z();


                DQDGamma_.boundaryField()[patchI][faceI] = //CQ.value();
                    tensor
                    (
                        DSDEpsilon.xx(), 0, 0,
                        0, CQ.value().yy(), 0,
                        0, 0, CQ.value().zz()
                    );

                DMDK_.boundaryField()[patchI][faceI] = //CM.value();
                    tensor
                    (
                        CM.value().xx(), 0,               0,
                        0,               DSDEpsilon.yy(), DSDEpsilon.yz(),
                        0,               DSDEpsilon.zy(), DSDEpsilon.zz()
                    );
                
                DMDGamma_.boundaryField()[patchI][faceI] = //tensor::zero;
                    tensor
                    (
                        0,               0, 0,
                        DSDEpsilon.yx(), 0, 0,
                        DSDEpsilon.zx(), 0, 0
                    );
                
                DQDK_.boundaryField()[patchI][faceI] = //tensor::zero;
                    tensor
                    (
                        0, DSDEpsilon.xy(), DSDEpsilon.xz(),
                        0, 0, 0,
                        0, 0, 0
                    );
                
                // if (patchI == 0 && faceI == 0)
                // {
                //     Info << "_Ftrial: " << Ftrial << endl;
                //     Info << "_F: " << F << endl;
                //     Info << "_NP: " << NP << endl;
                //     Info << "_N: " << N << endl;
                //     Info << "_DCP: " << DCP << endl;
                //     Info << "_DEP: " << DEP << endl;
                //     Info << "_dEP: " << dEP << endl;
                //     Info << "_DEpsilonP: " << DEpsilonP << endl;
                //     Info << "_dEpsilonP: " << dEpsilonP << endl;
                //     Info << "_______________DSDEpsilon: " << DSDEpsilon << endl;
                //     Info << "_CQ: " << CQ.value() << endl;
                //     Info << "_CM: " << CM.value() << endl;
                //     Info << "_DQDGamma_[faceI]: "
                //          << DQDGamma_.boundaryField()[patchI][faceI]
                //          << endl;
                //     Info << "_DMDK_[faceI]: "
                //          << DMDK_.boundaryField()[patchI][faceI] << endl;
                //     Info << "_My: " << MI[faceI].y() << endl;
                //     Info << "_Mz: " << MI[faceI].z() << endl;
                //     Info << "_M: " << M << endl;
                // }

                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const processorFvPatch& procPatch =
                        refCast<const processorFvPatch>(mesh.boundary()[patchI]);
                    if (procPatch.master())
                    {
                        nYieldingFaces_++;
                    }
                }
                else
                {
                    activeYield_.boundaryField()[patchI][faceI] = 1;
                    nYieldingFaces_++;
                }
                activeYield_[faceCells[faceI]] = 1;
            }
            else
            {
                // DGammaP_.boundaryField()[patchI][faceI] = vector::zero;
                // DKP_.boundaryField()[patchI][faceI] = vector::zero;
                    
                DGammaP_.boundaryField()[patchI][faceI] =
                (
                    oldGammaPI[faceI]
                  - GammaPI[faceI]
                );
                DKP_.boundaryField()[patchI][faceI] =
                (
                    oldKPI[faceI]
                  - KPI[faceI]
                );

                DQDGamma_.boundaryField()[patchI][faceI] = CQ.value();
                DMDK_.boundaryField()[patchI][faceI] = CM.value();
                
                DMDGamma_.boundaryField()[patchI][faceI] = tensor::zero;
                DQDK_.boundaryField()[patchI][faceI] = tensor::zero;
            }
        }
    }

    reduce(nYieldingFaces_, sumOp<label>());

    // Info << "nPlasticFaces: " << nPlasticFaces << endl;
}

void axialForce3dBending::updateYieldStress()
{
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

void axialForce3dBending::writeInfo() const
{
    Info << "Number of yielding faces: "
         << nYieldingFaces_ << endl;
}

} // end of namespace
// ************************************************************************* //
