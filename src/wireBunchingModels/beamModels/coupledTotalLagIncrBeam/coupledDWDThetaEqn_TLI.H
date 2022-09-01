{
    // Info << "Using coupled solver for DW and DTheta" << endl;

    scalar DThetaResidual = GREAT;
    scalar DWResidual = GREAT;
  
    DTheta_.storePrevIter();
    DW_.storePrevIter();

    // Initialise the block system
    fvBlockMatrix<vector6> DWDThetaEqn(DWDTheta_);

    // Grab block diagonal and set it to zero
    tensor6Field& d = DWDThetaEqn.diag().asSquare();
    d = tensor6::zero;

    // Grab off-diagonal and set it to zero
    tensor6Field& u = DWDThetaEqn.upper().asSquare();
    u = tensor6::zero;

    // Grab off-diagonal and set it to zero
    tensor6Field& l = DWDThetaEqn.lower().asSquare();
    l = tensor6::zero;

    // Grap source and set it to zero
    vector6Field& source = DWDThetaEqn.source();
    source = vector6::zero;

    scalarField deltaf = 1.0/mesh().deltaCoeffs().internalField();
    const scalarField& wf = mesh().weights().internalField();
    const vectorField& iI = i_.internalField();
    
    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    #include "updateCoefficients_TLI.H"

    const tensorField& CQDWI = CQDW_.internalField();
    const tensorField& CQDThetaI = CQDTheta_.internalField();
    const tensorField& CMDThetaI = CMDTheta_.internalField();

    const vectorField& explicitQI = explicitQ_.internalField();
    const vectorField& explicitMI = explicitM_.internalField();

    // surfaceVectorField Wf = fvc::interpolate(W_);

    // Internal faces
    forAll(u, faceI)
    {
        ////// W equation

        // W part (Laplacian)
        u[faceI](0,0) += CQDWI[faceI].xx()/deltaf[faceI];
        u[faceI](0,1) += CQDWI[faceI].xy()/deltaf[faceI];
        u[faceI](0,2) += CQDWI[faceI].xz()/deltaf[faceI];

        u[faceI](1,0) += CQDWI[faceI].yx()/deltaf[faceI];
        u[faceI](1,1) += CQDWI[faceI].yy()/deltaf[faceI];
        u[faceI](1,2) += CQDWI[faceI].yz()/deltaf[faceI];

        u[faceI](2,0) += CQDWI[faceI].zx()/deltaf[faceI];
        u[faceI](2,1) += CQDWI[faceI].zy()/deltaf[faceI];
        u[faceI](2,2) += CQDWI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](0,0) += -CQDWI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](0,1) += -CQDWI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](0,2) += -CQDWI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](1,0) += -CQDWI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](1,1) += -CQDWI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](1,2) += -CQDWI[faceI].yz()/deltaf[faceI];

        d[own[faceI]](2,0) += -CQDWI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](2,1) += -CQDWI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](2,2) += -CQDWI[faceI].zz()/deltaf[faceI];
        
        l[faceI](0,0) += CQDWI[faceI].xx()/deltaf[faceI];
        l[faceI](0,1) += CQDWI[faceI].xy()/deltaf[faceI];
        l[faceI](0,2) += CQDWI[faceI].xz()/deltaf[faceI];
        
        l[faceI](1,0) += CQDWI[faceI].yx()/deltaf[faceI];
        l[faceI](1,1) += CQDWI[faceI].yy()/deltaf[faceI];
        l[faceI](1,2) += CQDWI[faceI].yz()/deltaf[faceI];

        l[faceI](2,0) += CQDWI[faceI].zx()/deltaf[faceI];
        l[faceI](2,1) += CQDWI[faceI].zy()/deltaf[faceI];
        l[faceI](2,2) += CQDWI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](0,0) += -CQDWI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](0,1) += -CQDWI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](0,2) += -CQDWI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](1,0) += -CQDWI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](1,1) += -CQDWI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](1,2) += -CQDWI[faceI].yz()/deltaf[faceI];
        
        d[nei[faceI]](2,0) += -CQDWI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](2,1) += -CQDWI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](2,2) += -CQDWI[faceI].zz()/deltaf[faceI];

        // Theta part
        u[faceI](0,3) += (1-wf[faceI])*CQDThetaI[faceI].xx();
        u[faceI](0,4) += (1-wf[faceI])*CQDThetaI[faceI].xy();
        u[faceI](0,5) += (1-wf[faceI])*CQDThetaI[faceI].xz();

        u[faceI](1,3) += (1-wf[faceI])*CQDThetaI[faceI].yx();
        u[faceI](1,4) += (1-wf[faceI])*CQDThetaI[faceI].yy();
        u[faceI](1,5) += (1-wf[faceI])*CQDThetaI[faceI].yz();
        
        u[faceI](2,3) += (1-wf[faceI])*CQDThetaI[faceI].zx();
        u[faceI](2,4) += (1-wf[faceI])*CQDThetaI[faceI].zy();
        u[faceI](2,5) += (1-wf[faceI])*CQDThetaI[faceI].zz();
        
        d[own[faceI]](0,3) += wf[faceI]*CQDThetaI[faceI].xx();
        d[own[faceI]](0,4) += wf[faceI]*CQDThetaI[faceI].xy();
        d[own[faceI]](0,5) += wf[faceI]*CQDThetaI[faceI].xz();
        
        d[own[faceI]](1,3) += wf[faceI]*CQDThetaI[faceI].yx();
        d[own[faceI]](1,4) += wf[faceI]*CQDThetaI[faceI].yy();
        d[own[faceI]](1,5) += wf[faceI]*CQDThetaI[faceI].yz();
        
        d[own[faceI]](2,3) += wf[faceI]*CQDThetaI[faceI].zx();
        d[own[faceI]](2,4) += wf[faceI]*CQDThetaI[faceI].zy();
        d[own[faceI]](2,5) += wf[faceI]*CQDThetaI[faceI].zz();

        source[own[faceI]](0) -= explicitQI[faceI].x();
        source[own[faceI]](1) -= explicitQI[faceI].y();
        source[own[faceI]](2) -= explicitQI[faceI].z();

        l[faceI](0,3) += -wf[faceI]*CQDThetaI[faceI].xx();
        l[faceI](0,4) += -wf[faceI]*CQDThetaI[faceI].xy();
        l[faceI](0,5) += -wf[faceI]*CQDThetaI[faceI].xz();

        l[faceI](1,3) += -wf[faceI]*CQDThetaI[faceI].yx();
        l[faceI](1,4) += -wf[faceI]*CQDThetaI[faceI].yy();
        l[faceI](1,5) += -wf[faceI]*CQDThetaI[faceI].yz();

        l[faceI](2,3) += -wf[faceI]*CQDThetaI[faceI].zx();
        l[faceI](2,4) += -wf[faceI]*CQDThetaI[faceI].zy();
        l[faceI](2,5) += -wf[faceI]*CQDThetaI[faceI].zz();

        d[nei[faceI]](0,3) += -(1-wf[faceI])*CQDThetaI[faceI].xx();
        d[nei[faceI]](0,4) += -(1-wf[faceI])*CQDThetaI[faceI].xy();
        d[nei[faceI]](0,5) += -(1-wf[faceI])*CQDThetaI[faceI].xz();

        d[nei[faceI]](1,3) += -(1-wf[faceI])*CQDThetaI[faceI].yx();
        d[nei[faceI]](1,4) += -(1-wf[faceI])*CQDThetaI[faceI].yy();
        d[nei[faceI]](1,5) += -(1-wf[faceI])*CQDThetaI[faceI].yz();

        d[nei[faceI]](2,3) += -(1-wf[faceI])*CQDThetaI[faceI].zx();
        d[nei[faceI]](2,4) += -(1-wf[faceI])*CQDThetaI[faceI].zy();
        d[nei[faceI]](2,5) += -(1-wf[faceI])*CQDThetaI[faceI].zz();

        source[nei[faceI]](0) -= -explicitQI[faceI].x();
        source[nei[faceI]](1) -= -explicitQI[faceI].y();
        source[nei[faceI]](2) -= -explicitQI[faceI].z();

        
        ////// Theta equation

        // Laplacian part

        u[faceI](3,3) += CMDThetaI[faceI].xx()/deltaf[faceI];
        u[faceI](3,4) += CMDThetaI[faceI].xy()/deltaf[faceI];
        u[faceI](3,5) += CMDThetaI[faceI].xz()/deltaf[faceI];
        
        u[faceI](4,3) += CMDThetaI[faceI].yx()/deltaf[faceI];
        u[faceI](4,4) += CMDThetaI[faceI].yy()/deltaf[faceI];
        u[faceI](4,5) += CMDThetaI[faceI].yz()/deltaf[faceI];

        u[faceI](5,3) += CMDThetaI[faceI].zx()/deltaf[faceI];
        u[faceI](5,4) += CMDThetaI[faceI].zy()/deltaf[faceI];
        u[faceI](5,5) += CMDThetaI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](3,3) += -CMDThetaI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](3,4) += -CMDThetaI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](3,5) += -CMDThetaI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](4,3) += -CMDThetaI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](4,4) += -CMDThetaI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](4,5) += -CMDThetaI[faceI].yz()/deltaf[faceI];
        
        d[own[faceI]](5,3) += -CMDThetaI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](5,4) += -CMDThetaI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](5,5) += -CMDThetaI[faceI].zz()/deltaf[faceI];

        l[faceI](3,3) += CMDThetaI[faceI].xx()/deltaf[faceI];
        l[faceI](3,4) += CMDThetaI[faceI].xy()/deltaf[faceI];
        l[faceI](3,5) += CMDThetaI[faceI].xz()/deltaf[faceI];
        
        l[faceI](4,3) += CMDThetaI[faceI].yx()/deltaf[faceI];
        l[faceI](4,4) += CMDThetaI[faceI].yy()/deltaf[faceI];
        l[faceI](4,5) += CMDThetaI[faceI].yz()/deltaf[faceI];

        l[faceI](5,3) += CMDThetaI[faceI].zx()/deltaf[faceI];
        l[faceI](5,4) += CMDThetaI[faceI].zy()/deltaf[faceI];
        l[faceI](5,5) += CMDThetaI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](3,3) += -CMDThetaI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](3,4) += -CMDThetaI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](3,5) += -CMDThetaI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](4,3) += -CMDThetaI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](4,4) += -CMDThetaI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](4,5) += -CMDThetaI[faceI].yz()/deltaf[faceI];

        d[nei[faceI]](5,3) += -CMDThetaI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](5,4) += -CMDThetaI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](5,5) += -CMDThetaI[faceI].zz()/deltaf[faceI];

        //---- (dr x  Q) part

        vector curT = (Lambda_[faceI] & iI[faceI]);
        scalar curDelta = stretchRatio_[faceI]*deltaf[faceI];
        tensor DRP = spinTensor(wf[faceI]*curDelta*curT);
        tensor DRN = -spinTensor((1.0-wf[faceI])*curDelta*curT);
 
        // tensor DRP = spinTensor(Cf[faceI]-C[own[faceI]]);
        // tensor DRN = spinTensor(Cf[faceI]-C[nei[faceI]]);

        tensor DRPCQ = (DRP & CQDWI[faceI]);
        tensor DRNCQ = (DRN & CQDWI[faceI]);

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

        tensor DRPCQCT = ((DRP & CQDThetaI[faceI]));
        tensor DRNCQCT = ((DRN & CQDThetaI[faceI]));

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
    DW_.boundaryField().updateCoeffs();
    DTheta_.boundaryField().updateCoeffs();
    forAll (DW_.boundaryField(), patchI)
    {
        const tensorField& pCQDW = CQDW_.boundaryField()[patchI];
        const tensorField& pCQDTheta = CQDTheta_.boundaryField()[patchI];
        const tensorField& pCMDTheta = CMDTheta_.boundaryField()[patchI];

        const vectorField& pExplicitQ = explicitQ_.boundaryField()[patchI];
        const vectorField& pExplicitM = explicitM_.boundaryField()[patchI];
        
        const fvPatch& patch = mesh().boundary()[patchI];
        const labelList& fc = patch.faceCells();

        ////// DW equation

        if
        (
            isA<forceBeamDisplacementFvPatchVectorField>
            (
                DW_.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementFvPatchVectorField& pDW =
                refCast<forceBeamDisplacementFvPatchVectorField>
                (
                    DW_.boundaryField()[patchI]
                );

            // Source contribution
            forAll (pDW, faceI)
            {
                source[fc[faceI]](0) -= pDW.force()[faceI].x();
                source[fc[faceI]](1) -= pDW.force()[faceI].y();
                source[fc[faceI]](2) -= pDW.force()[faceI].z();
            }
        }
        else if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                DW_.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementFvPatchVectorField& pDW =
                refCast<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    DW_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];

            vectorField tang =
            (
                Lambda_.boundaryField()[patchI]
              & i_.boundaryField()[patchI]
            );

            tensorField pCQDWt = pCQDW - ((tang*tang) & pCQDW);
            tensorField pCQDThetat = pCQDTheta - ((tang*tang) & pCQDTheta);

            // Diag contribution
            forAll (pDW, faceI)
            {
                d[fc[faceI]](0,0) += -pCQDWt[faceI].xx()/pDelta[faceI];
                d[fc[faceI]](0,1) += -pCQDWt[faceI].xy()/pDelta[faceI];
                d[fc[faceI]](0,2) += -pCQDWt[faceI].xz()/pDelta[faceI];

                d[fc[faceI]](1,0) += -pCQDWt[faceI].yx()/pDelta[faceI];
                d[fc[faceI]](1,1) += -pCQDWt[faceI].yy()/pDelta[faceI];
                d[fc[faceI]](1,2) += -pCQDWt[faceI].yz()/pDelta[faceI];

                d[fc[faceI]](2,0) += -pCQDWt[faceI].zx()/pDelta[faceI];
                d[fc[faceI]](2,1) += -pCQDWt[faceI].zy()/pDelta[faceI];
                d[fc[faceI]](2,2) += -pCQDWt[faceI].zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pDW, faceI)
            {
                vector DWContrib =
                    (pCQDWt[faceI] & pDW.refDisp()[faceI])/pDelta[faceI];

                source[fc[faceI]](0) -= DWContrib.x();
                source[fc[faceI]](1) -= DWContrib.y();
                source[fc[faceI]](2) -= DWContrib.z();
            }

            if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pDTheta, faceI)
                {
                    vector thetaContrib =
                        (pCQDThetat[faceI] & pDTheta[faceI]);

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMDTheta/pDelta);

                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQDThetat[faceI]
                        & (invCM[faceI] & (pCMDTheta[faceI]/pDelta[faceI]))
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
                forAll (pDTheta, faceI)
                {
                    tensor CqCt = pCQDThetat[faceI];
                    vector thetaContrib =
                        (
                            CqCt
                          & (
                                invCM[faceI]
                              & (
                                    pDTheta.moment()[faceI]
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
            vectorField pQa = tang*pDW.axialForce();
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
                DW_.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pDW =
                refCast<fixedValueFvPatchVectorField>
                (
                    DW_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];

            // Diag contribution
            forAll (pDW, faceI)
            {
                d[fc[faceI]](0,0) += -pCQDW[faceI].xx()/pDelta[faceI];
                d[fc[faceI]](0,1) += -pCQDW[faceI].xy()/pDelta[faceI];
                d[fc[faceI]](0,2) += -pCQDW[faceI].xz()/pDelta[faceI];

                d[fc[faceI]](1,0) += -pCQDW[faceI].yx()/pDelta[faceI];
                d[fc[faceI]](1,1) += -pCQDW[faceI].yy()/pDelta[faceI];
                d[fc[faceI]](1,2) += -pCQDW[faceI].yz()/pDelta[faceI];

                d[fc[faceI]](2,0) += -pCQDW[faceI].zx()/pDelta[faceI];
                d[fc[faceI]](2,1) += -pCQDW[faceI].zy()/pDelta[faceI];
                d[fc[faceI]](2,2) += -pCQDW[faceI].zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pDW, faceI)
            {
                vector WContrib = (pCQDW[faceI] & pDW[faceI])/pDelta[faceI];
              
                source[fc[faceI]](0) -= WContrib.x();
                source[fc[faceI]](1) -= WContrib.y();
                source[fc[faceI]](2) -= WContrib.z();
            }

            if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pDTheta, faceI)
                {
                    vector thetaContrib =
                        (pCQDTheta[faceI] & pDTheta[faceI]);

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMDTheta/pDelta);
                
                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQDTheta[faceI]
                      & (invCM[faceI] & (pCMDTheta[faceI]/pDelta[faceI]))
                    );
                    
                    // tensor CqCt = 
                    //     (t.value() & pN[faceI])*pCQDTheta[faceI];

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
                forAll (pDTheta, faceI)
                {
                    tensor CqCt = pCQDTheta[faceI];
                    vector thetaContrib =
                        (
                            CqCt
                          & (
                                invCM[faceI]
                              & (
                                    pDTheta.moment()[faceI]
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
                DTheta_.boundaryField()[patchI]
            )
        )
        {
            const momentBeamRotationFvPatchVectorField& pDTheta =
                refCast<momentBeamRotationFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                );

            // Source contribution
            forAll (pDTheta, faceI)
            {
                source[fc[faceI]](3) -= pDTheta.moment()[faceI].x();
                source[fc[faceI]](4) -= pDTheta.moment()[faceI].y();
                source[fc[faceI]](5) -= pDTheta.moment()[faceI].z();
            }
        }
        else if
        (
            isA<fixedValueFvPatchVectorField>
            (
                DTheta_.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pDTheta =
                refCast<fixedValueFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];
            
            // Diag contribution from laplacian
            forAll (pDTheta, faceI)
            {
                d[fc[faceI]](3,3) += -pCMDTheta[faceI].xx()/pDelta[faceI];
                d[fc[faceI]](3,4) += -pCMDTheta[faceI].xy()/pDelta[faceI];
                d[fc[faceI]](3,5) += -pCMDTheta[faceI].xz()/pDelta[faceI];

                d[fc[faceI]](4,3) += -pCMDTheta[faceI].yx()/pDelta[faceI];
                d[fc[faceI]](4,4) += -pCMDTheta[faceI].yy()/pDelta[faceI];
                d[fc[faceI]](4,5) += -pCMDTheta[faceI].yz()/pDelta[faceI];

                d[fc[faceI]](5,3) += -pCMDTheta[faceI].zx()/pDelta[faceI];
                d[fc[faceI]](5,4) += -pCMDTheta[faceI].zy()/pDelta[faceI];
                d[fc[faceI]](5,5) += -pCMDTheta[faceI].zz()/pDelta[faceI];
            }

            // Source contribution
            forAll (pDTheta, faceI)
            {
                source[fc[faceI]](3) -=
                    (
                        pCMDTheta[faceI].xx()*pDTheta[faceI].x()
                      + pCMDTheta[faceI].xy()*pDTheta[faceI].y()
                      + pCMDTheta[faceI].xz()*pDTheta[faceI].z()
                    )
                   /pDelta[faceI];
                
                source[fc[faceI]](4) -=
                    (
                        pCMDTheta[faceI].yx()*pDTheta[faceI].x()
                      + pCMDTheta[faceI].yy()*pDTheta[faceI].y()
                      + pCMDTheta[faceI].yz()*pDTheta[faceI].z()
                    )
                   /pDelta[faceI];
                
                source[fc[faceI]](5) -=
                    (
                        pCMDTheta[faceI].zx()*pDTheta[faceI].x()
                      + pCMDTheta[faceI].zy()*pDTheta[faceI].y()
                      + pCMDTheta[faceI].zz()*pDTheta[faceI].z()
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
                DW_.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementFvPatchVectorField& pDW =
                refCast<forceBeamDisplacementFvPatchVectorField>
                (
                    DW_.boundaryField()[patchI]
                );

            // Source contribution
            forAll (pDW, faceI)
            {
                vector MFb = (DRb[faceI] & pDW.force()[faceI]);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }
        }
        else if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                DW_.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementFvPatchVectorField& pDW =
                refCast<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    DW_.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh().deltaCoeffs().boundaryField()[patchI];

            vectorField tang =
            (
                Lambda_.boundaryField()[patchI]
              & i_.boundaryField()[patchI]
            );

            tensorField pCQDWt = pCQDW - ((tang*tang) & pCQDW);
            tensorField pCQDThetat = pCQDTheta - ((tang*tang) & pCQDTheta);

            // Diag contribution
            forAll (pDW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQDWt[faceI]);

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
            forAll (pDW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQDWt[faceI]);
                vector curSource = (DRbCQ & pDW.refDisp()[faceI])/pDelta[faceI];

                source[fc[faceI]](3) -= curSource.x();
                source[fc[faceI]](4) -= curSource.y();
                source[fc[faceI]](5) -= curSource.z();
            }

            if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pDTheta, faceI)
                {
                    tensor DRbCQCT = (DRb[faceI] & pCQDThetat[faceI]);

                    vector curSource =
                        (DRbCQCT & pDTheta[faceI]);
                
                    source[fc[faceI]](3) -= curSource.x();
                    source[fc[faceI]](4) -= curSource.y();
                    source[fc[faceI]](5) -= curSource.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMDTheta/pDelta);

                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQDThetat[faceI]
                      & (invCM[faceI] & (pCMDTheta[faceI]/pDelta[faceI]))
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
                forAll (pDTheta, faceI)
                {
                    vector thetaContrib =
                    (
                        pCQDThetat[faceI]
                      & (
                            invCM[faceI]
                          & (
                                pDTheta.moment()[faceI]
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
            forAll (pDW, faceI)
            {
                vector pQa = tang[faceI]*pDW.axialForce()[faceI];
                vector MFb = (DRb[faceI] & pQa);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }

            // Add explicit force contribution
            forAll (pDW, faceI)
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
                DW_.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pDW =
                refCast<fixedValueFvPatchVectorField>
                (
                    DW_.boundaryField()[patchI]
                );

            // Diag contribution
            forAll (pDW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQDW[faceI]);

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
            forAll (pDW, faceI)
            {
                tensor DRbCQ = (DRb[faceI] & pCQDW[faceI]);
                vector curSource = (DRbCQ & pDW[faceI])/pDelta[faceI];

                source[fc[faceI]](3) -= curSource.x();
                source[fc[faceI]](4) -= curSource.y();
                source[fc[faceI]](5) -= curSource.z();
            }

            if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pDTheta, faceI)
                {
                    tensor DRbCQCT = (DRb[faceI] & pCQDTheta[faceI]);
                    vector curSource = (DRbCQCT & pDTheta[faceI]);
                
                    source[fc[faceI]](3) -= curSource.x();
                    source[fc[faceI]](4) -= curSource.y();
                    source[fc[faceI]](5) -= curSource.z();
                }
            }
            else if
            (
                isA<momentBeamRotationFvPatchVectorField>
                (
                    DTheta_.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta_.boundaryField()[patchI]
                    );

                scalarField pDelta =
                    1.0/mesh().deltaCoeffs().boundaryField()[patchI];

                tensorField invCM = inv(pCMDTheta/pDelta);

                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor DRbCQCT =
                    (
                        DRb[faceI]
                      & (
                            pCQDTheta[faceI]
                          & (
                                invCM[faceI]
                              & (pCMDTheta[faceI]/pDelta[faceI])
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
                forAll (pDTheta, faceI)
                {
                    tensor DRbCQCT = (DRb[faceI] & pCQDTheta[faceI]);

                    vector MFb =
                        (
                            DRbCQCT
                          & (
                                invCM[faceI]
                              & (
                                    pDTheta.moment()[faceI]
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


    // Block coupled solver call
    vector6 eqnRes = DWDThetaEqn.solve().initialResidual();

    currentResidual = mag(eqnRes);

    // if (iOuterCorr() == 0)
    // {
    //     initialResidual = currentResidual;
    // }
    
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

    DWDThetaEqn.retrieveSolution(3, DTheta_.internalField());
    DTheta_.correctBoundaryConditions();

    DWDThetaEqn.retrieveSolution(0, DW_.internalField());
    DW_.correctBoundaryConditions();

    Q_.storePrevIter();
    Q_ = (CQDW_ & fvc::snGrad(DW_))
      + (CQDTheta_ & fvc::interpolate(DTheta_))
      + explicitQ_;
    Q_.relax();

    M_.storePrevIter();
    M_ = (CMDTheta_ & fvc::snGrad(DTheta_))
      + (CMDTheta2_ & fvc::interpolate(DTheta_))
      + explicitM_;
    M_.relax();

    DW_.relax();
    DTheta_.relax();

    // Insert under-relaxed DTheta field into DWDTheta block field
    if (false)
    {
        const direction nCmpts = pTraits<vector>::nComponents;
        direction localDir = 3;

        Field<vector6>& psiIn = DWDTheta_.internalField();

        for (direction cmptI = 0; cmptI < nCmpts; cmptI++)
        {
            scalarField xSingleCurr(DTheta_.component(cmptI));

            forAll (xSingleCurr, cellI)
            {
                psiIn[cellI](localDir) = xSingleCurr[cellI];
            }

            localDir++;
        }
    }

    // Calculate axial stretch ratio
    if (true)
    {
        surfaceVectorField Gamma =
        (
            Lambda_.T()
          & (
                fvc::snGrad(W_.oldTime() + DW_)
              + ((refLambda_ - Lambda_) & i_)
            )
        );

        surfaceScalarField GammaA = Gamma.component(0);
        GammaA.boundaryField()[startPatchIndex()] =
           -GammaA.boundaryField()[startPatchIndex()];

        stretchRatio_.storePrevIter();
        stretchRatio_ = mag(1 + GammaA);
        stretchRatio_.relax();
    }

    // Calculate DTheta residual
    if (false)
    {
        scalar denom =
            gMax
            (
                mag
                (
                    DTheta_.internalField()
                )
            );

        if (denom < SMALL)
        {
            denom = max(gMax(mag(DTheta_.internalField())), SMALL);
        }

        DThetaResidual =
            gMax
            (
                mag
                (
                    DTheta_.internalField()
                  - DTheta_.prevIter().internalField()
                )
            )/denom;
    }

    // Calculate DW residual
    if (false)
    {
        scalar denom =
            gMax
            (
                mag
                (
                    DW_.internalField()
                )
            );

        if (denom < SMALL)
        {
            denom = max(gMax(mag(DW_.internalField())), SMALL);
        }

        DWResidual =
            gMax
            (
                mag
                (
                    DW_.internalField()
                  - DW_.prevIter().internalField()
                )
            )/denom;
    }
    max(DWResidual, DThetaResidual);
 
    // currentResidual = max(DWResidual, DThetaResidual);

    if (iOuterCorr() == 0)
    // if ( (iCorr == 0) || (currentResidual > initialResidual) )
    {
        initialResidual = currentResidual;
    }

    // Info << currentResidual << endl;
}
