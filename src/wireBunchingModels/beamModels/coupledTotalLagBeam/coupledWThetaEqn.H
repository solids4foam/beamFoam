{
    // Info << "Using coupled solver for DW and DTheta" << endl;

    // scalar ThetaResidual = GREAT;
    // scalar WResidual = GREAT;
  
    Theta_.storePrevIter();
    // W_.storePrevIter();

    // Initialise the block system
    fvBlockMatrix<vector6> WThetaEqn(WTheta_);

    // Grab block diagonal and set it to zero
    tensor6Field& d = WThetaEqn.diag().asSquare();
    d = tensor6::zero;

    // Grab off-diagonal and set it to zero
    tensor6Field& u = WThetaEqn.upper().asSquare();
    u = tensor6::zero;

    // Grab off-diagonal and set it to zero
    tensor6Field& l = WThetaEqn.lower().asSquare();
    l = tensor6::zero;

    // Grap source and set it to zero
    vector6Field& source = WThetaEqn.source();
    source = vector6::zero;

    scalarField deltaf = 1.0/mesh().deltaCoeffs().internalField();
    const scalarField& wf = mesh().weights().internalField();
    const vectorField& iI = i_.internalField();
    
    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    #include "updateCoefficients_TL.H"

    const tensorField& CQWI = CQW_.internalField();
    const tensorField& CQThetaI = CQTheta_.internalField();
    const tensorField& CMThetaI = CMTheta_.internalField();

    const vectorField& explicitQI = explicitQ_.internalField();
    const vectorField& explicitMI = explicitM_.internalField();

    surfaceVectorField Wf = fvc::interpolate(W_);

    // Internal faces
    forAll(u, faceI)
    {
        ////// W equation

        // W part (Laplacian)
        u[faceI](0,0) += CQWI[faceI].xx()/deltaf[faceI];
        u[faceI](0,1) += CQWI[faceI].xy()/deltaf[faceI];
        u[faceI](0,2) += CQWI[faceI].xz()/deltaf[faceI];

        u[faceI](1,0) += CQWI[faceI].yx()/deltaf[faceI];
        u[faceI](1,1) += CQWI[faceI].yy()/deltaf[faceI];
        u[faceI](1,2) += CQWI[faceI].yz()/deltaf[faceI];

        u[faceI](2,0) += CQWI[faceI].zx()/deltaf[faceI];
        u[faceI](2,1) += CQWI[faceI].zy()/deltaf[faceI];
        u[faceI](2,2) += CQWI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](1,2) += -CQWI[faceI].yz()/deltaf[faceI];

        d[own[faceI]](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](2,2) += -CQWI[faceI].zz()/deltaf[faceI];
        
        l[faceI](0,0) += CQWI[faceI].xx()/deltaf[faceI];
        l[faceI](0,1) += CQWI[faceI].xy()/deltaf[faceI];
        l[faceI](0,2) += CQWI[faceI].xz()/deltaf[faceI];
        
        l[faceI](1,0) += CQWI[faceI].yx()/deltaf[faceI];
        l[faceI](1,1) += CQWI[faceI].yy()/deltaf[faceI];
        l[faceI](1,2) += CQWI[faceI].yz()/deltaf[faceI];

        l[faceI](2,0) += CQWI[faceI].zx()/deltaf[faceI];
        l[faceI](2,1) += CQWI[faceI].zy()/deltaf[faceI];
        l[faceI](2,2) += CQWI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](1,2) += -CQWI[faceI].yz()/deltaf[faceI];
        
        d[nei[faceI]](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](2,2) += -CQWI[faceI].zz()/deltaf[faceI];

        // Theta part
        u[faceI](0,3) += (1-wf[faceI])*CQThetaI[faceI].xx();
        u[faceI](0,4) += (1-wf[faceI])*CQThetaI[faceI].xy();
        u[faceI](0,5) += (1-wf[faceI])*CQThetaI[faceI].xz();

        u[faceI](1,3) += (1-wf[faceI])*CQThetaI[faceI].yx();
        u[faceI](1,4) += (1-wf[faceI])*CQThetaI[faceI].yy();
        u[faceI](1,5) += (1-wf[faceI])*CQThetaI[faceI].yz();
        
        u[faceI](2,3) += (1-wf[faceI])*CQThetaI[faceI].zx();
        u[faceI](2,4) += (1-wf[faceI])*CQThetaI[faceI].zy();
        u[faceI](2,5) += (1-wf[faceI])*CQThetaI[faceI].zz();
        
        d[own[faceI]](0,3) += wf[faceI]*CQThetaI[faceI].xx();
        d[own[faceI]](0,4) += wf[faceI]*CQThetaI[faceI].xy();
        d[own[faceI]](0,5) += wf[faceI]*CQThetaI[faceI].xz();
        
        d[own[faceI]](1,3) += wf[faceI]*CQThetaI[faceI].yx();
        d[own[faceI]](1,4) += wf[faceI]*CQThetaI[faceI].yy();
        d[own[faceI]](1,5) += wf[faceI]*CQThetaI[faceI].yz();
        
        d[own[faceI]](2,3) += wf[faceI]*CQThetaI[faceI].zx();
        d[own[faceI]](2,4) += wf[faceI]*CQThetaI[faceI].zy();
        d[own[faceI]](2,5) += wf[faceI]*CQThetaI[faceI].zz();

        source[own[faceI]](0) -= explicitQI[faceI].x();
        source[own[faceI]](1) -= explicitQI[faceI].y();
        source[own[faceI]](2) -= explicitQI[faceI].z();

        l[faceI](0,3) += -wf[faceI]*CQThetaI[faceI].xx();
        l[faceI](0,4) += -wf[faceI]*CQThetaI[faceI].xy();
        l[faceI](0,5) += -wf[faceI]*CQThetaI[faceI].xz();

        l[faceI](1,3) += -wf[faceI]*CQThetaI[faceI].yx();
        l[faceI](1,4) += -wf[faceI]*CQThetaI[faceI].yy();
        l[faceI](1,5) += -wf[faceI]*CQThetaI[faceI].yz();

        l[faceI](2,3) += -wf[faceI]*CQThetaI[faceI].zx();
        l[faceI](2,4) += -wf[faceI]*CQThetaI[faceI].zy();
        l[faceI](2,5) += -wf[faceI]*CQThetaI[faceI].zz();

        d[nei[faceI]](0,3) += -(1-wf[faceI])*CQThetaI[faceI].xx();
        d[nei[faceI]](0,4) += -(1-wf[faceI])*CQThetaI[faceI].xy();
        d[nei[faceI]](0,5) += -(1-wf[faceI])*CQThetaI[faceI].xz();

        d[nei[faceI]](1,3) += -(1-wf[faceI])*CQThetaI[faceI].yx();
        d[nei[faceI]](1,4) += -(1-wf[faceI])*CQThetaI[faceI].yy();
        d[nei[faceI]](1,5) += -(1-wf[faceI])*CQThetaI[faceI].yz();

        d[nei[faceI]](2,3) += -(1-wf[faceI])*CQThetaI[faceI].zx();
        d[nei[faceI]](2,4) += -(1-wf[faceI])*CQThetaI[faceI].zy();
        d[nei[faceI]](2,5) += -(1-wf[faceI])*CQThetaI[faceI].zz();

        source[nei[faceI]](0) -= -explicitQI[faceI].x();
        source[nei[faceI]](1) -= -explicitQI[faceI].y();
        source[nei[faceI]](2) -= -explicitQI[faceI].z();

        
        ////// Theta equation

        // Laplacian part

        u[faceI](3,3) += CMThetaI[faceI].xx()/deltaf[faceI];
        u[faceI](3,4) += CMThetaI[faceI].xy()/deltaf[faceI];
        u[faceI](3,5) += CMThetaI[faceI].xz()/deltaf[faceI];
        
        u[faceI](4,3) += CMThetaI[faceI].yx()/deltaf[faceI];
        u[faceI](4,4) += CMThetaI[faceI].yy()/deltaf[faceI];
        u[faceI](4,5) += CMThetaI[faceI].yz()/deltaf[faceI];

        u[faceI](5,3) += CMThetaI[faceI].zx()/deltaf[faceI];
        u[faceI](5,4) += CMThetaI[faceI].zy()/deltaf[faceI];
        u[faceI](5,5) += CMThetaI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];
        
        d[own[faceI]](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];

        l[faceI](3,3) += CMThetaI[faceI].xx()/deltaf[faceI];
        l[faceI](3,4) += CMThetaI[faceI].xy()/deltaf[faceI];
        l[faceI](3,5) += CMThetaI[faceI].xz()/deltaf[faceI];
        
        l[faceI](4,3) += CMThetaI[faceI].yx()/deltaf[faceI];
        l[faceI](4,4) += CMThetaI[faceI].yy()/deltaf[faceI];
        l[faceI](4,5) += CMThetaI[faceI].yz()/deltaf[faceI];

        l[faceI](5,3) += CMThetaI[faceI].zx()/deltaf[faceI];
        l[faceI](5,4) += CMThetaI[faceI].zy()/deltaf[faceI];
        l[faceI](5,5) += CMThetaI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];

        d[nei[faceI]](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];

        //---- (dr x  Q) part

        vector curT = (Lambda_[faceI] & iI[faceI]);
        scalar curDelta = stretchRatio_[faceI]*deltaf[faceI];
        tensor DRP = spinTensor(wf[faceI]*curDelta*curT);
        tensor DRN = -spinTensor((1.0-wf[faceI])*curDelta*curT);
 
        // tensor DRP = spinTensor(Cf[faceI]-C[own[faceI]]);
        // tensor DRN = spinTensor(Cf[faceI]-C[nei[faceI]]);

        tensor DRPCQ = (DRP & CQWI[faceI]);
        tensor DRNCQ = (DRN & CQWI[faceI]);

        // W part
        u[faceI](3,0) += DRPCQ.xx()/deltaf[faceI];
        u[faceI](3,1) += DRPCQ.xy()/deltaf[faceI];
        u[faceI](3,2) += DRPCQ.xz()/deltaf[faceI];
        
        u[faceI](4,0) += DRPCQ.yx()/deltaf[faceI];
        u[faceI](4,1) += DRPCQ.yy()/deltaf[faceI];
        u[faceI](4,2) += DRPCQ.yz()/deltaf[faceI];
        
        u[faceI](5,0) += DRPCQ.zx()/deltaf[faceI];
        u[faceI](5,1) += DRPCQ.zy()/deltaf[faceI];
        u[faceI](5,2) += DRPCQ.zz()/deltaf[faceI];

        d[own[faceI]](3,0) -= DRPCQ.xx()/deltaf[faceI];
        d[own[faceI]](3,1) -= DRPCQ.xy()/deltaf[faceI];
        d[own[faceI]](3,2) -= DRPCQ.xz()/deltaf[faceI];

        d[own[faceI]](4,0) -= DRPCQ.yx()/deltaf[faceI];
        d[own[faceI]](4,1) -= DRPCQ.yy()/deltaf[faceI];
        d[own[faceI]](4,2) -= DRPCQ.yz()/deltaf[faceI];
        
        d[own[faceI]](5,0) -= DRPCQ.zx()/deltaf[faceI];
        d[own[faceI]](5,1) -= DRPCQ.zy()/deltaf[faceI];
        d[own[faceI]](5,2) -= DRPCQ.zz()/deltaf[faceI];

        l[faceI](3,0) += DRNCQ.xx()/deltaf[faceI];
        l[faceI](3,1) += DRNCQ.xy()/deltaf[faceI];
        l[faceI](3,2) += DRNCQ.xz()/deltaf[faceI];

        l[faceI](4,0) += DRNCQ.yx()/deltaf[faceI];
        l[faceI](4,1) += DRNCQ.yy()/deltaf[faceI];
        l[faceI](4,2) += DRNCQ.yz()/deltaf[faceI];

        l[faceI](5,0) += DRNCQ.zx()/deltaf[faceI];
        l[faceI](5,1) += DRNCQ.zy()/deltaf[faceI];
        l[faceI](5,2) += DRNCQ.zz()/deltaf[faceI];

        d[nei[faceI]](3,0) -= DRNCQ.xx()/deltaf[faceI];
        d[nei[faceI]](3,1) -= DRNCQ.xy()/deltaf[faceI];
        d[nei[faceI]](3,2) -= DRNCQ.xz()/deltaf[faceI];

        d[nei[faceI]](4,0) -= DRNCQ.yx()/deltaf[faceI];
        d[nei[faceI]](4,1) -= DRNCQ.yy()/deltaf[faceI];
        d[nei[faceI]](4,2) -= DRNCQ.yz()/deltaf[faceI];

        d[nei[faceI]](5,0) -= DRNCQ.zx()/deltaf[faceI];
        d[nei[faceI]](5,1) -= DRNCQ.zy()/deltaf[faceI];
        d[nei[faceI]](5,2) -= DRNCQ.zz()/deltaf[faceI];


        // Theta part

        tensor DRPCQCT = ((DRP & CQThetaI[faceI]));
        tensor DRNCQCT = ((DRN & CQThetaI[faceI]));

        u[faceI](3,3) += (1-wf[faceI])*DRPCQCT.xx();
        u[faceI](3,4) += (1-wf[faceI])*DRPCQCT.xy();
        u[faceI](3,5) += (1-wf[faceI])*DRPCQCT.xz();

        u[faceI](4,3) += (1-wf[faceI])*DRPCQCT.yx();
        u[faceI](4,4) += (1-wf[faceI])*DRPCQCT.yy();
        u[faceI](4,5) += (1-wf[faceI])*DRPCQCT.yz();

        u[faceI](5,3) += (1-wf[faceI])*DRPCQCT.zx();
        u[faceI](5,4) += (1-wf[faceI])*DRPCQCT.zy();
        u[faceI](5,5) += (1-wf[faceI])*DRPCQCT.zz();

        d[own[faceI]](3,3) += wf[faceI]*DRPCQCT.xx();
        d[own[faceI]](3,4) += wf[faceI]*DRPCQCT.xy();
        d[own[faceI]](3,5) += wf[faceI]*DRPCQCT.xz();

        d[own[faceI]](4,3) += wf[faceI]*DRPCQCT.yx();
        d[own[faceI]](4,4) += wf[faceI]*DRPCQCT.yy();
        d[own[faceI]](4,5) += wf[faceI]*DRPCQCT.yz();
        
        d[own[faceI]](5,3) += wf[faceI]*DRPCQCT.zx();
        d[own[faceI]](5,4) += wf[faceI]*DRPCQCT.zy();
        d[own[faceI]](5,5) += wf[faceI]*DRPCQCT.zz();
        
        l[faceI](3,3) -= wf[faceI]*DRNCQCT.xx();
        l[faceI](3,4) -= wf[faceI]*DRNCQCT.xy();
        l[faceI](3,5) -= wf[faceI]*DRNCQCT.xz();
        
        l[faceI](4,3) -= wf[faceI]*DRNCQCT.yx();
        l[faceI](4,4) -= wf[faceI]*DRNCQCT.yy();
        l[faceI](4,5) -= wf[faceI]*DRNCQCT.yz();
        
        l[faceI](5,3) -= wf[faceI]*DRNCQCT.zx();
        l[faceI](5,4) -= wf[faceI]*DRNCQCT.zy();
        l[faceI](5,5) -= wf[faceI]*DRNCQCT.zz();

        d[nei[faceI]](3,3) -= (1-wf[faceI])*DRNCQCT.xx();
        d[nei[faceI]](3,4) -= (1-wf[faceI])*DRNCQCT.xy();
        d[nei[faceI]](3,5) -= (1-wf[faceI])*DRNCQCT.xz();
        
        d[nei[faceI]](4,3) -= (1-wf[faceI])*DRNCQCT.yx();
        d[nei[faceI]](4,4) -= (1-wf[faceI])*DRNCQCT.yy();
        d[nei[faceI]](4,5) -= (1-wf[faceI])*DRNCQCT.yz();
        
        d[nei[faceI]](5,3) -= (1-wf[faceI])*DRNCQCT.zx();
        d[nei[faceI]](5,4) -= (1-wf[faceI])*DRNCQCT.zy();
        d[nei[faceI]](5,5) -= (1-wf[faceI])*DRNCQCT.zz();


        // Explicit force
        vector explicitQMP = (DRP & explicitQI[faceI]);
        source[own[faceI]](3) -= explicitQMP.x();
        source[own[faceI]](4) -= explicitQMP.y();
        source[own[faceI]](5) -= explicitQMP.z();

        vector explicitQMN = (DRN & explicitQI[faceI]);
        source[nei[faceI]](3) -= -explicitQMN.x();
        source[nei[faceI]](4) -= -explicitQMN.y();
        source[nei[faceI]](5) -= -explicitQMN.z();

        source[own[faceI]](3) -= explicitMI[faceI].x();
        source[own[faceI]](4) -= explicitMI[faceI].y();
        source[own[faceI]](5) -= explicitMI[faceI].z();
        
        source[nei[faceI]](3) -= -explicitMI[faceI].x();
        source[nei[faceI]](4) -= -explicitMI[faceI].y();
        source[nei[faceI]](5) -= -explicitMI[faceI].z();
    }

    // Boundary contributions
    W_.boundaryField().updateCoeffs();
    Theta_.boundaryField().updateCoeffs();
    forAll (W_.boundaryField(), patchI)
    {
        const tensorField& pCQW = CQW_.boundaryField()[patchI];
        const tensorField& pCQTheta = CQTheta_.boundaryField()[patchI];
        const tensorField& pCMTheta = CMTheta_.boundaryField()[patchI];

        const vectorField& pExplicitQ = explicitQ_.boundaryField()[patchI];
        const vectorField& pExplicitM = explicitM_.boundaryField()[patchI];

        const fvPatch& patch = mesh().boundary()[patchI];
        const labelList& fc = patch.faceCells();

        ////// W equation

        if
        (
            isA<forceBeamDisplacementFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementFvPatchVectorField& pW =
                refCast<forceBeamDisplacementFvPatchVectorField>
                (
                    W_.boundaryField()[patchI]
                );

            // Source contribution
            forAll (pW, faceI)
            {
                source[fc[faceI]](0) -= pW.force()[faceI].x();
                source[fc[faceI]](1) -= pW.force()[faceI].y();
                source[fc[faceI]](2) -= pW.force()[faceI].z();
            }
        }
        else if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementFvPatchVectorField& pW =
                refCast<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    W_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];

            vectorField tang =
            (
                Lambda_.boundaryField()[patchI]
              & i_.boundaryField()[patchI]
            );

            tensorField pCQWt = pCQW - ((tang*tang) & pCQW);
            tensorField pCQThetat = pCQTheta - ((tang*tang) & pCQTheta);

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
                    (pCQWt[faceI] & pW.refDisp()[faceI])/pDelta[faceI];

                source[fc[faceI]](0) -= WContrib.x();
                source[fc[faceI]](1) -= WContrib.y();
                source[fc[faceI]](2) -= WContrib.z();
            }

            if
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
                        Theta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector thetaContrib =
                        (pCQThetat[faceI] & pTheta[faceI]);

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        Theta_.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMTheta/pDelta);

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

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }

            // Add axial component of force
            vectorField pQa = tang*pW.axialForce();
            forAll (pQa, faceI)
            {
                source[fc[faceI]](0) -= pQa[faceI].x();
                source[fc[faceI]](1) -= pQa[faceI].y();
                source[fc[faceI]](2) -= pQa[faceI].z();
            }
            
            // Add explicit force contribution
            vectorField pExplicitQt = pExplicitQ - ((tang*tang) & pExplicitQ);
            forAll (pExplicitQ, faceI)
            {
                source[fc[faceI]](0) -= pExplicitQt[faceI].x();
                source[fc[faceI]](1) -= pExplicitQt[faceI].y();
                source[fc[faceI]](2) -= pExplicitQt[faceI].z();
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
                    W_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];

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
                        pW[faceI]
                      + vector(SMALL, SMALL, SMALL) // Not sure why
                    )
                )/pDelta[faceI];
              
                source[fc[faceI]](0) -= WContrib.x();
                source[fc[faceI]](1) -= WContrib.y();
                source[fc[faceI]](2) -= WContrib.z();
            }

            if
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
                        Theta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector thetaContrib =
                        (pCQTheta[faceI] & pTheta[faceI]);

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        Theta_.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMTheta/pDelta);
                
                // Diag contribution
                forAll (pTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQTheta[faceI]
                      & (invCM[faceI] & (pCMTheta[faceI]/pDelta[faceI]))
                    );
                    
                    // tensor CqCt = 
                    //     (t.value() & pN[faceI])*pCQTheta[faceI];

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

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }

            // Add explicit force conribution
            forAll (pExplicitQ, faceI)
            {
                source[fc[faceI]](0) -= pExplicitQ[faceI].x();
                source[fc[faceI]](1) -= pExplicitQ[faceI].y();
                source[fc[faceI]](2) -= pExplicitQ[faceI].z();
            }
        }

        ////// Theta equation
        
        if
        (
            isA<momentBeamRotationFvPatchVectorField>
            (
                Theta_.boundaryField()[patchI]
            )
        )
        {
            const momentBeamRotationFvPatchVectorField& pTheta =
                refCast<momentBeamRotationFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                );

            // Source contribution
            forAll (pTheta, faceI)
            {
                source[fc[faceI]](3) -= pTheta.moment()[faceI].x();
                source[fc[faceI]](4) -= pTheta.moment()[faceI].y();
                source[fc[faceI]](5) -= pTheta.moment()[faceI].z();
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
                    Theta_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];
            
            // Diag contribution from laplacian
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
                source[fc[faceI]](3) -=
                    (
                        pCMTheta[faceI].xx()*pTheta[faceI].x()
                      + pCMTheta[faceI].xy()*pTheta[faceI].y()
                      + pCMTheta[faceI].xz()*pTheta[faceI].z()
                    )
                   /pDelta[faceI];
                
                source[fc[faceI]](4) -=
                    (
                        pCMTheta[faceI].yx()*pTheta[faceI].x()
                      + pCMTheta[faceI].yy()*pTheta[faceI].y()
                      + pCMTheta[faceI].yz()*pTheta[faceI].z()
                    )
                   /pDelta[faceI];
                
                source[fc[faceI]](5) -=
                    (
                        pCMTheta[faceI].zx()*pTheta[faceI].x()
                      + pCMTheta[faceI].zy()*pTheta[faceI].y()
                      + pCMTheta[faceI].zz()*pTheta[faceI].z()
                    )
                   /pDelta[faceI];

                // Explicit moment
                source[fc[faceI]](3) -= pExplicitM[faceI].x();
                source[fc[faceI]](4) -= pExplicitM[faceI].y();
                source[fc[faceI]](5) -= pExplicitM[faceI].z();
            }
        }

        // dr x Q term

        scalarField pDelta =
            1.0/mesh().deltaCoeffs().boundaryField()[patchI];

        const vectorField& pi = i_.boundaryField()[patchI];

        vectorField t =
        (
            Lambda_.boundaryField()[patchI] & pi
        );

        vectorField curDelta =
            stretchRatio_.boundaryField()[patchI]*pDelta*t;

        tensorField DRb = spinTensor(curDelta);

        if
        (
            isA<forceBeamDisplacementFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementFvPatchVectorField& pW =
                refCast<forceBeamDisplacementFvPatchVectorField>
                (
                    W_.boundaryField()[patchI]
                );

            // Source contribution
            forAll (pW, faceI)
            {
                vector MFb = (DRb[faceI] & pW.force()[faceI]);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }
        }
        else if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                W_.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementFvPatchVectorField& pW =
                refCast<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    W_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];

            vectorField tang =
            (
                Lambda_.boundaryField()[patchI]
              & i_.boundaryField()[patchI]
            );

            tensorField pCQWt = pCQW - ((tang*tang) & pCQW);
            tensorField pCQThetat = pCQTheta - ((tang*tang) & pCQTheta);

            // Diag contribution
            forAll (pW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQWt[faceI]);

                d[fc[faceI]](3,0) -= DRbCQ.xx()/pDelta[faceI];
                d[fc[faceI]](3,1) -= DRbCQ.xy()/pDelta[faceI];
                d[fc[faceI]](3,2) -= DRbCQ.xz()/pDelta[faceI];

                d[fc[faceI]](4,0) -= DRbCQ.yx()/pDelta[faceI];
                d[fc[faceI]](4,1) -= DRbCQ.yy()/pDelta[faceI];
                d[fc[faceI]](4,2) -= DRbCQ.yz()/pDelta[faceI];

                d[fc[faceI]](5,0) -= DRbCQ.zx()/pDelta[faceI];
                d[fc[faceI]](5,1) -= DRbCQ.zy()/pDelta[faceI];
                d[fc[faceI]](5,2) -= DRbCQ.zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQWt[faceI]);
                vector curSource = (DRbCQ & pW.refDisp()[faceI])/pDelta[faceI];

                source[fc[faceI]](3) -= curSource.x();
                source[fc[faceI]](4) -= curSource.y();
                source[fc[faceI]](5) -= curSource.z();
            }

            if
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
                        Theta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pTheta, faceI)
                {
                    tensor DRbCQCT = (DRb[faceI] & pCQThetat[faceI]);

                    vector curSource =
                        (DRbCQCT & pTheta[faceI]);
                
                    source[fc[faceI]](3) -= curSource.x();
                    source[fc[faceI]](4) -= curSource.y();
                    source[fc[faceI]](5) -= curSource.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        Theta_.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMTheta/pDelta);

                // Diag contribution
                forAll (pTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQThetat[faceI]
                      & (invCM[faceI] & (pCMTheta[faceI]/pDelta[faceI]))
                    );
                    tensor DRbCQCT = (DRb[faceI] & CqCt);

                    d[fc[faceI]](3,3) += DRbCQCT.xx();
                    d[fc[faceI]](3,4) += DRbCQCT.xy();
                    d[fc[faceI]](3,5) += DRbCQCT.xz();

                    d[fc[faceI]](4,3) += DRbCQCT.yx();
                    d[fc[faceI]](4,4) += DRbCQCT.yy();
                    d[fc[faceI]](4,5) += DRbCQCT.yz();

                    d[fc[faceI]](5,3) += DRbCQCT.zx();
                    d[fc[faceI]](5,4) += DRbCQCT.zy();
                    d[fc[faceI]](5,5) += DRbCQCT.zz();
                }

                // Source contribution
                forAll (pTheta, faceI)
                {
                    vector thetaContrib =
                    (
                        pCQThetat[faceI]
                      & (
                            invCM[faceI]
                          & (
                                pTheta.moment()[faceI]
                              - pExplicitM[faceI]
                            )
                        )
                    );
                    
                    vector MFb = (DRb[faceI] & thetaContrib);

                    source[fc[faceI]](3) -= MFb.x();
                    source[fc[faceI]](4) -= MFb.y();
                    source[fc[faceI]](5) -= MFb.z();
                }
            }
            
            // Add axial force
            forAll (pW, faceI)
            {
                vector pQa = tang[faceI]*pW.axialForce()[faceI];
                vector MFb = (DRb[faceI] & pQa);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }

            // Add explicit force contribution
            forAll (pW, faceI)
            {
                vector pQt = pExplicitQ[faceI];
                pQt -= ((tang[faceI]*tang[faceI]) & pQt);

                vector MFb = (DRb[faceI] & pQt);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
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
                    W_.boundaryField()[patchI]
                );

            // Diag contribution
            forAll (pW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQW[faceI]);

                d[fc[faceI]](3,0) -= DRbCQ.xx()/pDelta[faceI];
                d[fc[faceI]](3,1) -= DRbCQ.xy()/pDelta[faceI];
                d[fc[faceI]](3,2) -= DRbCQ.xz()/pDelta[faceI];

                d[fc[faceI]](4,0) -= DRbCQ.yx()/pDelta[faceI];
                d[fc[faceI]](4,1) -= DRbCQ.yy()/pDelta[faceI];
                d[fc[faceI]](4,2) -= DRbCQ.yz()/pDelta[faceI];

                d[fc[faceI]](5,0) -= DRbCQ.zx()/pDelta[faceI];
                d[fc[faceI]](5,1) -= DRbCQ.zy()/pDelta[faceI];
                d[fc[faceI]](5,2) -= DRbCQ.zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQW[faceI]);
                vector curSource =
                (
                    DRbCQ
                  & (
                        pW[faceI]
                      // + vector(SMALL, SMALL, SMALL)
                    )
                )/pDelta[faceI];

                source[fc[faceI]](3) -= curSource.x();
                source[fc[faceI]](4) -= curSource.y();
                source[fc[faceI]](5) -= curSource.z();
            }

            if
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
                        Theta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pTheta, faceI)
                {
                    tensor DRbCQCT = (DRb[faceI] & pCQTheta[faceI]);
                    vector curSource = (DRbCQCT & pTheta[faceI]);
                
                    source[fc[faceI]](3) -= curSource.x();
                    source[fc[faceI]](4) -= curSource.y();
                    source[fc[faceI]](5) -= curSource.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    Theta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        Theta_.boundaryField()[patchI]
                    );

                scalarField pDelta =
                    1.0/mesh().deltaCoeffs().boundaryField()[patchI];

                tensorField invCM = inv(pCMTheta/pDelta);

                // Diag contribution
                forAll (pTheta, faceI)
                {
                    tensor DRbCQCT =
                    (
                        DRb[faceI]
                      & (
                            pCQTheta[faceI]
                          & (
                                invCM[faceI]
                              & (pCMTheta[faceI]/pDelta[faceI])
                            )
                        )
                    );

                    d[fc[faceI]](3,3) += DRbCQCT.xx();
                    d[fc[faceI]](3,4) += DRbCQCT.xy();
                    d[fc[faceI]](3,5) += DRbCQCT.xz();

                    d[fc[faceI]](4,3) += DRbCQCT.yx();
                    d[fc[faceI]](4,4) += DRbCQCT.yy();
                    d[fc[faceI]](4,5) += DRbCQCT.yz();

                    d[fc[faceI]](5,3) += DRbCQCT.zx();
                    d[fc[faceI]](5,4) += DRbCQCT.zy();
                    d[fc[faceI]](5,5) += DRbCQCT.zz();
                }

                // Source contribution
                forAll (pTheta, faceI)
                {
                    tensor DRbCQCT = (DRb[faceI] & pCQTheta[faceI]);

                    vector MFb =
                        (
                            DRbCQCT
                          & (
                                invCM[faceI]
                              & (
                                    pTheta.moment()[faceI]
                                  - pExplicitM[faceI]
                                  // - pOldM[faceI]
                                  // - pDM[faceI]
                                )
                            )
                        );

                    source[fc[faceI]](3) -= MFb.x();
                    source[fc[faceI]](4) -= MFb.y();
                    source[fc[faceI]](5) -= MFb.z();
                }
            }

            // Add explicit force contribution
            forAll (pExplicitQ, faceI)
            {
                vector MFb = (DRb[faceI] & pExplicitQ[faceI]);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }
        }
    }

    // Add distributed force
    forAll(source, cellI)
    {
        source[cellI](0) -= q()[cellI].x()*L()[cellI];
        source[cellI](1) -= q()[cellI].y()*L()[cellI];
        source[cellI](2) -= q()[cellI].z()*L()[cellI];
    }

    // Add distributed moment
    forAll(source, cellI)
    {
        source[cellI](3) -= m()[cellI].x()*L()[cellI];
        source[cellI](4) -= m()[cellI].y()*L()[cellI];
        source[cellI](5) -= m()[cellI].z()*L()[cellI];
    }
    
    // Block coupled solver call
    vector6 eqnRes = WThetaEqn.solve().initialResidual();

    currentResidual = mag(eqnRes);

    if (iOuterCorr() == 0)
    {
        initialResidual = currentResidual;
    }
    
    // maxResidual =
    //     max
    //     (
    //         mag
    //         (
    //             DWDTheta.internalField()
    //           - DWDTheta.prevIter().internalField()
    //         )
    //        /max(mag(DWDTheta.internalField()))
    //     );

    // Info << "Max residual: " << maxResidual << endl;

    WThetaEqn.retrieveSolution(3, Theta_.internalField());
    Theta_.correctBoundaryConditions();

    WThetaEqn.retrieveSolution(0, W_.internalField());
    W_.correctBoundaryConditions();

    epsilon_ = fvc::snGrad(W_) + ((I - Lambda_) & i_);

    Q_.storePrevIter();
    Q_ = (CQW_ & fvc::snGrad(W_))
      + (CQTheta_ & fvc::interpolate(Theta_))
      + explicitQ_;
    Q_.relax();

    M_.storePrevIter();
    M_ = (CMTheta_ & fvc::snGrad(Theta_)) + explicitM_;
    M_.relax();

    // W_.relax();
    Theta_.relax();

    // Calculate axial stretch ratio
    if (true)
    {
        surfaceVectorField Epsilon =
        (
            Lambda_.T()
          & (
                fvc::snGrad(W_)
              + ((refLambda_ - Lambda_) & i_)
            )
        );

        surfaceScalarField EpsilonA = Epsilon.component(0);
        EpsilonA.boundaryField()[startPatchIndex()] =
           -EpsilonA.boundaryField()[startPatchIndex()];

        stretchRatio_.storePrevIter();
        stretchRatio_ = mag(1 + EpsilonA);
        stretchRatio_.relax();
    }

    // // Calculate DTheta residual
    // {
    //     scalar denom =
    //         gMax
    //         (
    //             mag
    //             (
    //                 Theta_.internalField()
    //               - Theta_.oldTime().internalField()
    //             )
    //         );

    //     if (denom < SMALL)
    //     {
    //         denom = max(gMax(mag(Theta_.internalField())), SMALL);
    //     }

    //     ThetaResidual =
    //         gMax
    //         (
    //             mag
    //             (
    //                 Theta_.internalField()
    //               - Theta_.prevIter().internalField()
    //             )
    //         )/denom;
    // }

    // // Calculate DW residual
    // {
    //     scalar denom =
    //         gMax
    //         (
    //             mag
    //             (
    //                 W_.internalField()
    //               - W_.oldTime().internalField()
    //             )
    //         );

    //     if (denom < SMALL)
    //     {
    //         denom = max(gMax(mag(W_.internalField())), SMALL);
    //     }

    //     WResidual =
    //         gMax
    //         (
    //             mag
    //             (
    //                 W_.internalField()
    //               - W_.prevIter().internalField()
    //             )
    //         )/denom;
    // }

    // currentResidual = max(WResidual, ThetaResidual);

    // if (iCorr == 0)
    // {
    //     initialResidual = currentResidual;
    // }
}
