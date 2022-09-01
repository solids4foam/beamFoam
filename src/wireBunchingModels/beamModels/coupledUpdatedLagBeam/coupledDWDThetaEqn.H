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
    const vectorField& dRuDsI = dRuDs_.internalField();

    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();
    // const vectorField& Cf = mesh().Cf().internalField();
    // const vectorField& C = mesh().C().internalField();

    #include "updateCoefficients.H"

    const tensorField& CQDWI = CQDW_.internalField();
    const tensorField& CQDThetaI = CQDTheta_.internalField();
    const tensorField& CMDThetaI = CMDTheta_.internalField();
    const tensorField& CMDTheta2I = CMDTheta2_.internalField();

    // const vectorField& dQI = dQ_.internalField();
    // const vectorField& QI = Q_.internalField();
    const vectorField& oldQI = Q_.oldTime().internalField();

    const vectorField& oldMI = M_.oldTime().internalField();
    // const vectorField& dMI = dM.internalField();

    const vectorField& DWI = DW_.internalField();
    surfaceVectorField DWf = fvc::interpolate(DW_);
    // const vectorField& DWfI = DWf.internalField();

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

        // // Non-linear correction
        // source[own[faceI]](0) -= dQI[faceI].x();
        // source[own[faceI]](1) -= dQI[faceI].y();
        // source[own[faceI]](2) -= dQI[faceI].z();

        vector oldQ = oldQI[faceI];
        // vector oldQ = (DLambdaI[faceI] & oldQI[faceI]);
        source[own[faceI]](0) -= oldQ.x();
        source[own[faceI]](1) -= oldQ.y();
        source[own[faceI]](2) -= oldQ.z();

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

        // // Non-linear correction
        // source[nei[faceI]](0) -= -dQI[faceI].x();
        // source[nei[faceI]](1) -= -dQI[faceI].y();
        // source[nei[faceI]](2) -= -dQI[faceI].z();

        source[nei[faceI]](0) -= -oldQ.x();
        source[nei[faceI]](1) -= -oldQ.y();
        source[nei[faceI]](2) -= -oldQ.z();
        
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

        vector curT =
        (
            DLambda_[faceI]
          & dRuDsI[faceI]
        );
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


        // Implicit old moment term
        u[faceI](3,3) += (1-wf[faceI])*CMDTheta2I[faceI].xx();
        u[faceI](3,4) += (1-wf[faceI])*CMDTheta2I[faceI].xy();
        u[faceI](3,5) += (1-wf[faceI])*CMDTheta2I[faceI].xz();

        u[faceI](4,3) += (1-wf[faceI])*CMDTheta2I[faceI].yx();
        u[faceI](4,4) += (1-wf[faceI])*CMDTheta2I[faceI].yy();
        u[faceI](4,5) += (1-wf[faceI])*CMDTheta2I[faceI].yz();
        
        u[faceI](5,3) += (1-wf[faceI])*CMDTheta2I[faceI].zx();
        u[faceI](5,4) += (1-wf[faceI])*CMDTheta2I[faceI].zy();
        u[faceI](5,5) += (1-wf[faceI])*CMDTheta2I[faceI].zz();
        
        d[own[faceI]](3,3) += wf[faceI]*CMDTheta2I[faceI].xx();
        d[own[faceI]](3,4) += wf[faceI]*CMDTheta2I[faceI].xy();
        d[own[faceI]](3,5) += wf[faceI]*CMDTheta2I[faceI].xz();

        d[own[faceI]](4,3) += wf[faceI]*CMDTheta2I[faceI].yx();
        d[own[faceI]](4,4) += wf[faceI]*CMDTheta2I[faceI].yy();
        d[own[faceI]](4,5) += wf[faceI]*CMDTheta2I[faceI].yz();

        d[own[faceI]](5,3) += wf[faceI]*CMDTheta2I[faceI].zx();
        d[own[faceI]](5,4) += wf[faceI]*CMDTheta2I[faceI].zy();
        d[own[faceI]](5,5) += wf[faceI]*CMDTheta2I[faceI].zz();

        l[faceI](3,3) += -wf[faceI]*CMDTheta2I[faceI].xx();
        l[faceI](3,4) += -wf[faceI]*CMDTheta2I[faceI].xy();
        l[faceI](3,5) += -wf[faceI]*CMDTheta2I[faceI].xz();

        l[faceI](4,3) += -wf[faceI]*CMDTheta2I[faceI].yx();
        l[faceI](4,4) += -wf[faceI]*CMDTheta2I[faceI].yy();
        l[faceI](4,5) += -wf[faceI]*CMDTheta2I[faceI].yz();

        l[faceI](5,3) += -wf[faceI]*CMDTheta2I[faceI].zx();
        l[faceI](5,4) += -wf[faceI]*CMDTheta2I[faceI].zy();
        l[faceI](5,5) += -wf[faceI]*CMDTheta2I[faceI].zz();
        
        d[nei[faceI]](3,3) += -(1-wf[faceI])*CMDTheta2I[faceI].xx();
        d[nei[faceI]](3,4) += -(1-wf[faceI])*CMDTheta2I[faceI].xy();
        d[nei[faceI]](3,5) += -(1-wf[faceI])*CMDTheta2I[faceI].xz();
        
        d[nei[faceI]](4,3) += -(1-wf[faceI])*CMDTheta2I[faceI].yx();
        d[nei[faceI]](4,4) += -(1-wf[faceI])*CMDTheta2I[faceI].yy();
        d[nei[faceI]](4,5) += -(1-wf[faceI])*CMDTheta2I[faceI].yz();
        
        d[nei[faceI]](5,3) += -(1-wf[faceI])*CMDTheta2I[faceI].zx();
        d[nei[faceI]](5,4) += -(1-wf[faceI])*CMDTheta2I[faceI].zy();
        d[nei[faceI]](5,5) += -(1-wf[faceI])*CMDTheta2I[faceI].zz();


        
        // // Nonlinear correction
        // vector MP = (DRP & dQI[faceI]);
        // source[own[faceI]](3) -= MP.x();
        // source[own[faceI]](4) -= MP.y();
        // source[own[faceI]](5) -= MP.z();
        
        // vector MN = (DRN & dQI[faceI]);
        // source[nei[faceI]](3) -= -MN.x();
        // source[nei[faceI]](4) -= -MN.y();
        // source[nei[faceI]](5) -= -MN.z();


        // // ZT
        // // Nonlinear correction (dW)
        // vector dWxdQP = ((DWfI[faceI] - DWI[own[faceI]]) ^ QI[faceI]);
        // source[own[faceI]](3) -= dWxdQP.x();
        // source[own[faceI]](4) -= dWxdQP.y();
        // source[own[faceI]](5) -= dWxdQP.z();

        // vector dWxdQN = ((DWfI[faceI] - DWI[nei[faceI]]) ^ QI[faceI]);
        // source[nei[faceI]](3) -= -dWxdQN.x();
        // source[nei[faceI]](4) -= -dWxdQN.y();
        // source[nei[faceI]](5) -= -dWxdQN.z();


        // Old time force
        vector oldQMP = (DRP & oldQI[faceI]);
        // vector oldQMP = (DRP & (DLambdaI[faceI] & oldQI[faceI]));
        source[own[faceI]](3) -= oldQMP.x();
        source[own[faceI]](4) -= oldQMP.y();
        source[own[faceI]](5) -= oldQMP.z();

        vector oldQMN = (DRN & oldQI[faceI]);
        // vector oldQMN = (DRN & (DLambdaI[faceI] & oldQI[faceI]));
        source[nei[faceI]](3) -= -oldQMN.x();
        source[nei[faceI]](4) -= -oldQMN.y();
        source[nei[faceI]](5) -= -oldQMN.z();

        // Old time moment and non-linear correction
        vector oldM = oldMI[faceI];
        // vector oldM = oldMI[faceI] + dMI[faceI];
        // vector oldM = (DLambdaI[faceI] & oldMI[faceI]);
        source[own[faceI]](3) -= oldM.x();
        source[own[faceI]](4) -= oldM.y();
        source[own[faceI]](5) -= oldM.z();
        
        source[nei[faceI]](3) -= -oldM.x();
        source[nei[faceI]](4) -= -oldM.y();
        source[nei[faceI]](5) -= -oldM.z();
    }

    if (false) //ZT
    {
        forAll (DW_.boundaryField(), patchI)
        {
            const fvPatch& patch = mesh().boundary()[patchI];
            const labelList& fc = patch.faceCells();
            const vectorField& pQ = Q_.boundaryField()[patchI];
            const vectorField& pDWf = DWf.boundaryField()[patchI];

            // const labelListList& cCells = mesh.cellCells();
            
            if
            (
                isA<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    DW_.boundaryField()[patchI]
                )
            )
            {
                // Add non-linear correcton of force (dW)
                forAll (pQ, faceI)
                {
                    label P = fc[faceI];
                    vector DWP = DWI[P];
                    // label PP = cCells[P][0];
                    // vector DWPP = DWI[PP];
              
                    vector DWb = DWP;
                    // vector DWb = 3*DWP/2 - DWPP/2;
              
                    vector MFb = ((DWb - DWP) ^ pQ[faceI]);
            
                    source[fc[faceI]](3) -= MFb.x();
                    source[fc[faceI]](4) -= MFb.y();
                    source[fc[faceI]](5) -= MFb.z();
                }
            }
            else
            {
                // Add non-linear correcton of force (dW)
                forAll (pQ, faceI)
                {
                    vector MFb =
                        ((pDWf[faceI] - DWI[fc[faceI]]) ^ pQ[faceI]);

                    source[fc[faceI]](3) -= MFb.x();
                    source[fc[faceI]](4) -= MFb.y();
                    source[fc[faceI]](5) -= MFb.z();
                }
            }
        }
    }

    // Boundary contributions
    DW_.boundaryField().updateCoeffs();
    DTheta_.boundaryField().updateCoeffs();
    
    forAll (DW_.boundaryField(), patchI)
    {
        const tensorField& pCQDW = CQDW_.boundaryField()[patchI];
        const tensorField& pCQDTheta = CQDTheta_.boundaryField()[patchI];
        const tensorField& pCMDTheta = CMDTheta_.boundaryField()[patchI];
        const tensorField& pCMDTheta2 = CMDTheta2_.boundaryField()[patchI];

        // const vectorField& pCf =
        //     mesh().Cf().boundaryField()[patchI];

        const vectorField& pOldQ = Q_.oldTime().boundaryField()[patchI];
        //   + QTheta.boundaryField()[patchI];

        // const vectorField& pDM = dM.boundaryField()[patchI];
        const vectorField& pOldM = M_.oldTime().boundaryField()[patchI];
        // const vectorField& pInitM = initM.boundaryField()[patchI];

        const fvPatch& patch = mesh().boundary()[patchI];
        const labelList& fc = patch.faceCells();

        ////// W equation

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
                DLambda_.boundaryField()[patchI]
              & dRuDs_.boundaryField()[patchI]
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

                tensorField invCM = inv(pCMDTheta/pDelta + pCMDTheta2);

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
                                  - pOldM[faceI]
                                  // - pOldM[faceI] - pDM[faceI]
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
            
            // // Add non-linear correcton of force
            // vectorField pDQt = pDQ - ((tang*tang) & pDQ);
            // forAll (pDQ, faceI)
            // {
            //     source[fc[faceI]](0) -= pDQt[faceI].x();
            //     source[fc[faceI]](1) -= pDQt[faceI].y();
            //     source[fc[faceI]](2) -= pDQt[faceI].z();
            // }

            // // Add non-linear correcton of force (old time)
            // vectorField pOldQt =
            //     pOldQ
            //   - ((tang*tang) & pOldQ);
            // forAll (pOldQ, faceI)
            // {
            //     // vector oldQ = pOldQt[faceI];
            //     source[fc[faceI]](0) -= pOldQt[faceI].x();
            //     source[fc[faceI]](1) -= pOldQt[faceI].y();
            //     source[fc[faceI]](2) -= pOldQt[faceI].z();
            // }

            // Add non-linear correcton of force (old time)
            vectorField pOldQt =
                pOldQ
              - ((tang*tang) & pOldQ);
            forAll (pOldQ, faceI)
            {
                // vector oldQ = pOldQt[faceI];
                source[fc[faceI]](0) -= pOldQt[faceI].x();
                source[fc[faceI]](1) -= pOldQt[faceI].y();
                source[fc[faceI]](2) -= pOldQt[faceI].z();
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
                vector DWContrib = (pCQDW[faceI] & pDW[faceI])/pDelta[faceI];
              
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

                tensorField invCM = inv(pCMDTheta/pDelta + pCMDTheta2);
                
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
                                  - pOldM[faceI]
                                  // - pOldM[faceI] - pDM[faceI]
                                )
                            )
                        );

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }

            // // Add non-linear correcton of force
            // forAll (pDQ, faceI)
            // {
            //     source[fc[faceI]](0) -= pDQ[faceI].x();
            //     source[fc[faceI]](1) -= pDQ[faceI].y();
            //     source[fc[faceI]](2) -= pDQ[faceI].z();
            // }
            
            // Add non-linear correcton of force (old time)
            forAll (pOldQ, faceI)
            {
                vector oldQ = pOldQ[faceI];
                source[fc[faceI]](0) -= oldQ.x();
                source[fc[faceI]](1) -= oldQ.y();
                source[fc[faceI]](2) -= oldQ.z();
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

                // Implicit old time moment term contribution
                source[fc[faceI]](3) -=
                    (
                        pCMDTheta2[faceI].xx()*pDTheta[faceI].x()
                      + pCMDTheta2[faceI].xy()*pDTheta[faceI].y()
                      + pCMDTheta2[faceI].xz()*pDTheta[faceI].z()
                    );
                
                source[fc[faceI]](4) -=
                    (
                        pCMDTheta2[faceI].yx()*pDTheta[faceI].x()
                      + pCMDTheta2[faceI].yy()*pDTheta[faceI].y()
                      + pCMDTheta2[faceI].yz()*pDTheta[faceI].z()
                    );
                
                source[fc[faceI]](5) -=
                    (
                        pCMDTheta2[faceI].zx()*pDTheta[faceI].x()
                      + pCMDTheta2[faceI].zy()*pDTheta[faceI].y()
                      + pCMDTheta2[faceI].zz()*pDTheta[faceI].z()
                    );

                // Old time moment
                vector oldM = pOldM[faceI];
                // vector oldM = pOldM[faceI] + pDM[faceI];
                // vector oldM = (pDLambda[faceI] & pOldM[faceI]);
                source[fc[faceI]](3) -= oldM.x();
                source[fc[faceI]](4) -= oldM.y();
                source[fc[faceI]](5) -= oldM.z();
            }
        }

        // dr x Q term

        scalarField pDelta =
            1.0/mesh().deltaCoeffs().boundaryField()[patchI];

        const vectorField& pDRuDs =
            dRuDs_.boundaryField()[patchI];

        vectorField t =
        (
            DLambda_.boundaryField()[patchI]
          & pDRuDs
        );

        vectorField curDelta = //(pCf - vectorField(C, fc));
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

            // scalarField pDelta =
            //     1.0/mesh().deltaCoeffs().boundaryField()[patchI];
            
            // const vectorField& pDRuDs =
            //     dRuDs_.boundaryField()[patchI];

            // Source contribution
            forAll (pDW, faceI)
            {
                // vector t =
                // (
                //     DLambda_.boundaryField()[patchI][faceI]
                //   & pDRuDs[faceI]
                // );
                // // vector t = pDRuDs[faceI];
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];                
                // tensor DRb = spinTensor(curDelta*t);
                // tensor DRb =
                //     spinTensor
                //     (
                //         pCf[faceI]
                //       - C[fc[faceI]]
                //     );
              
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
                DLambda_.boundaryField()[patchI]
              & dRuDs_.boundaryField()[patchI]
            );

            tensorField pCQDWt = pCQDW - ((tang*tang) & pCQDW);
            tensorField pCQDThetat = pCQDTheta - ((tang*tang) & pCQDTheta);

            // Diag contribution
            forAll (pDW, faceI)
            {
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];
                // tensor DRb = spinTensor(curDelta*tang[faceI]);
                // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
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
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];                
                // tensor DRb = spinTensor(curDelta*tang[faceI]);
                // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                tensor DRbCQ = (DRb[faceI] & pCQDWt[faceI]);

                vector curSource =
                    (DRbCQ & pDW.refDisp()[faceI])/pDelta[faceI];

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

                // vectorField tang =
                // (
                //     DLambda_.boundaryField()[patchI]
                //   & dRuDs_.boundaryField()[patchI]
                // );
            
                // vectorField pDThetat = pDTheta;
                // pDThetat -= tang*(tang & pDTheta);

                // Source contribution
                forAll (pDTheta, faceI)
                {
                    // scalar curDelta =
                    //     stretchRatio_.boundaryField()[patchI][faceI]
                    //    *pDelta[faceI];
                    // tensor DRb = spinTensor(curDelta*tang[faceI]);
                    // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
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

                // vectorField tang =
                // (
                //     DLambda_.boundaryField()[patchI]
                //   & dRuDs_.boundaryField()[patchI]
                // );
                
                tensorField invCM = inv(pCMDTheta/pDelta+pCMDTheta2);

                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    // scalar curDelta =
                    //     stretchRatio_.boundaryField()[patchI][faceI]
                    //    *pDelta[faceI];
                    // tensor DRb = spinTensor(curDelta*tang[faceI]);
                    // tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );
                    tensor CqCt =
                    (
                        pCQDThetat[faceI]
                      & (invCM[faceI] & (pCMDTheta[faceI]/pDelta[faceI]))
                    );
                    // CqCt -= (CqCt & (tang[faceI]*tang[faceI]));
                    tensor DRbCQCT = (DRb[faceI] & CqCt);
                    // (
                    //     DRb
                    //   & (
                    //         pCQDTheta[faceI]
                    //       & (
                    //             invCM[faceI]
                    //           & (pCMDTheta[faceI]/pDelta[faceI])
                    //         )
                    //     )
                    // );

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
                    // scalar curDelta =
                    //     stretchRatio_.boundaryField()[patchI][faceI]
                    //    *pDelta[faceI];
                    // tensor DRb = spinTensor(curDelta*tang[faceI]);
                    // tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );

                    // tensor CqCt = pCQDTheta[faceI];
                    vector thetaContrib =
                    (
                        // CqCt
                        pCQDThetat[faceI]
                      & (
                            invCM[faceI]
                          & (
                                pDTheta.moment()[faceI]
                              - pOldM[faceI]
                              // - pOldM[faceI] - pDM[faceI]
                            )
                        )
                    );
                    // thetaContrib -=
                    //     (thetaContrib & (tang[faceI]*tang[faceI]));
                    
                    // tensor DRbCQCT = (DRb & pCQDTheta[faceI]);
                    
                    vector MFb = (DRb[faceI] & thetaContrib);
                        // (
                        //     DRbCQCT
                        //   & (
                        //         invCM[faceI]
                        //       & (
                        //             pDTheta.moment()[faceI]
                        //           - pOldM[faceI]
                        //           - pDM[faceI]
                        //         )
                        //     )
                        // );

                    source[fc[faceI]](3) -= MFb.x();
                    source[fc[faceI]](4) -= MFb.y();
                    source[fc[faceI]](5) -= MFb.z();
                }
            }
            
            // Add axial force
            forAll (pDW, faceI)
            {
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];
                // tensor DRb = spinTensor(curDelta*tang[faceI]);
                // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                vector pQa = tang[faceI]*pDW.axialForce()[faceI];
                vector MFb = (DRb[faceI] & pQa);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }

            // Add non-linear correcton of force
            forAll (pDW, faceI)
            {
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];
                // tensor DRb = spinTensor(curDelta*tang[faceI]);
                // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);

                vector pQt = pOldQ[faceI];
                // vector pQt = pDQ[faceI] + pOldQ[faceI];
                pQt -= ((tang[faceI]*tang[faceI]) & pQt);

                vector MFb = (DRb[faceI] & pQt);
                    // (DRb & (pDQ[faceI] + pOldQ[faceI]));

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

            // scalarField pDelta =
            //     1.0/mesh().deltaCoeffs().boundaryField()[patchI];

            // const vectorField& pDRuDs = dRuDs_.boundaryField()[patchI];

            // Diag contribution
            forAll (pDW, faceI)
            {
                // vector t =
                // (
                //     DLambda_.boundaryField()[patchI][faceI]
                //   & pDRuDs[faceI]
                // );
                // vector t = pDRuDs[faceI];
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];                
                // tensor DRb = spinTensor(curDelta*t);
                // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
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
            //     vector t =
            //     (
            //         DLambda_.boundaryField()[patchI][faceI]
            //       & pDRuDs[faceI]
            //     );
                // vector t = pDRuDs[faceI];              
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];                
                // tensor DRb = spinTensor(curDelta*t);
                // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
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
                    // vector t =
                    // (
                    //     DLambda_.boundaryField()[patchI][faceI]
                    //   & pDRuDs[faceI]
                    // );
                    // scalar curDelta =
                    //     stretchRatio_.boundaryField()[patchI][faceI]
                    //    *pDelta[faceI];                
                    // tensor DRb = spinTensor(curDelta*t);
                    // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
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

                // const vectorField& pDRuDs = dRuDs_.boundaryField()[patchI];
            
                tensorField invCM = inv(pCMDTheta/pDelta+pCMDTheta2);

                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    // vector t =
                    // (
                    //     DLambda_.boundaryField()[patchI][faceI]
                    //   & pDRuDs[faceI]
                    // );
                    // vector t = pDRuDs[faceI];
                    // scalar curDelta =
                    //     stretchRatio_.boundaryField()[patchI][faceI]
                    //    *pDelta[faceI];                
                    // tensor DRb = spinTensor(curDelta*t);
                    // tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );
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
                    // vector t =
                    // (
                    //     DLambda_.boundaryField()[patchI][faceI]
                    //   & pDRuDs[faceI]
                    // );
                    // vector t = pDRuDs[faceI];
                    // scalar curDelta =
                    //     stretchRatio_.boundaryField()[patchI][faceI]
                    //    *pDelta[faceI];
                    // tensor DRb = spinTensor(curDelta*t);
                    // tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );
                    tensor DRbCQCT = (DRb[faceI] & pCQDTheta[faceI]);
                    
                    vector MFb =
                        (
                            DRbCQCT
                          & (
                                invCM[faceI]
                              & (
                                    pDTheta.moment()[faceI]
                                  - pOldM[faceI]
                                  // - pDM[faceI]
                                )
                            )
                        );

                    source[fc[faceI]](3) -= MFb.x();
                    source[fc[faceI]](4) -= MFb.y();
                    source[fc[faceI]](5) -= MFb.z();
                }
            }

            // Add non-linear correcton of force
            forAll (pOldQ, faceI)
            {
                // vector t =
                // (
                //     DLambda_.boundaryField()[patchI][faceI] & pDRuDs[faceI]
                // );
                // vector t = pDRuDs[faceI];
                // scalar curDelta =
                //     stretchRatio_.boundaryField()[patchI][faceI]
                //    *pDelta[faceI];                
                // tensor DRb = spinTensor(curDelta*t);
                // tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                vector MFb =
                    (DRb[faceI] & pOldQ[faceI]);
                    // (DRb & (pDQ[faceI] + pOldQ[faceI]));

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }        
        }

        // // Add non-linear correcton of moment (old time)
        // forAll (pOldM, faceI)
        // {
        //     vector oldM =
        //         (pDLambda[faceI] & pOldM[faceI]);
        //     // vector oldM =
        //     //     (t.value() & pN[faceI])
        //     //    *(pDLambda[faceI] & pOldM[faceI]);
        //     source[fc[faceI]](3) -= oldM.x();
        //     source[fc[faceI]](4) -= oldM.y();
        //     source[fc[faceI]](5) -= oldM.z();
        // }
    }

// Info << q() << endl;

    // Add distributed force
    forAll(source, cellI)
    {
        source[cellI](0) -= q()[cellI].x()*L()[cellI];
        source[cellI](1) -= q()[cellI].y()*L()[cellI];
        source[cellI](2) -= q()[cellI].z()*L()[cellI];
    }

    // scalar alpha = 0.9;
    // source += (d & DWDTheta)*(1-alpha)/alpha;
    // d = d/alpha;

    
    // Block coupled solver call
    vector6 eqnRes = DWDThetaEqn.solve().initialResidual();

    currentResidual = mag(eqnRes);

    if (iCorr == 0)
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

    DWDThetaEqn.retrieveSolution(3, DTheta_.internalField());
    DTheta_.correctBoundaryConditions();

    // DTheta.relax();

    // Info << "Calculate Q" << endl;
    DWDThetaEqn.retrieveSolution(0, DW_.internalField());
    // DW.boundaryField().updateCoeffs();
    DW_.correctBoundaryConditions();

    // Q = (CQDW & fvc::snGrad(DW)) + (CQDW & ((II - DLambda) & dRuDs))
    //   + (DLambda & Q.oldTime());

    // DW.relax();
    // DTheta.relax();

    epsilon_ =
        fvc::snGrad(DW_)
      + ((I - DLambda_) & dRuDs_)
      + (DLambda_ & epsilon_.oldTime());

    Q_.storePrevIter();
    Q_ = (CQDW_ & fvc::snGrad(DW_))
      + (CQDTheta_ & fvc::interpolate(DTheta_)) // + dQ
      + Q_.oldTime();
      // + (DLambda & Q.oldTime());
    Q_.relax();

    // Q_ = dimensionedVector("0", Q_.dimensions(), vector::zero);

    M_.storePrevIter();
    M_ = (CMDTheta_ & fvc::snGrad(DTheta_))
      + (CMDTheta2_ & fvc::interpolate(DTheta_)) // + dM
      + M_.oldTime();
      // + (DLambda & M.oldTime());
    M_.relax();

    // surfaceVectorField snGradDTheta = fvc::snGrad(DTheta_);
    // Info << snGradDTheta.boundaryField() << endl;
    // Info << snGradDTheta[0] << endl;
    
    // Info << DW.boundaryField() << endl;

    DW_.relax();
    DTheta_.relax();

    // Calculate axial stretch ratio
    if (true)
    {
        surfaceVectorField DEpsilon =
        (
            DLambda_.T()
          & (
                fvc::snGrad(DW_)
              + ((I - DLambda_) & dRuDs_)
            )
        );

        surfaceScalarField DEpsilonA = DEpsilon.component(0);
        DEpsilonA.boundaryField()[startPatchIndex()] =
           -DEpsilonA.boundaryField()[startPatchIndex()];

        stretchRatio_.storePrevIter();
        stretchRatio_ = mag(1 + DEpsilonA);
        stretchRatio_.relax();
    }

    // Calculate DTheta residual
    {
        scalar denom =
            gMax
            (
                mag
                (
                    DTheta_.internalField()
                  - DTheta_.oldTime().internalField()
                )
            );

        // if (denom < SMALL)
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
    {
        scalar denom =
            gMax
            (
                mag
                (
                    DW_.internalField()
                  - DW_.oldTime().internalField()
                )
            );

        // if (denom < SMALL)
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

    currentResidual = max(DWResidual, DThetaResidual);
    
    if (iCorr == 0)
    {
        initialResidual = currentResidual;
    }
    
    // DTheta_.internalField() = vector::zero;
}
