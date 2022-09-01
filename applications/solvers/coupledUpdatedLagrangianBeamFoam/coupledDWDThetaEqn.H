{
    Info << "Using coupled solver for DW and DTheta" << endl;

    DTheta.storePrevIter();
    DW.storePrevIter();
    // DWDTheta.storePrevIter();

    // Initialise the block system
    fvBlockMatrix<vector6> DWDThetaEqn(DWDTheta);

    // Grab block diagonal and set it to zero
    tensor6Field& d = DWDThetaEqn.diag().asSquare();
    d = tensor6::zero;

    // Info << DWDThetaEqn.diag().asSquare().size() << endl;

    // Grab off-diagonal and set it to zero
    tensor6Field& u = DWDThetaEqn.upper().asSquare();
    u = tensor6::zero;

    // Grab off-diagonal and set it to zero
    tensor6Field& l = DWDThetaEqn.lower().asSquare();
    l = tensor6::zero;

    // Grap source and set it to zero
    vector6Field& source = DWDThetaEqn.source();
    source = vector6::zero;

    scalarField deltaf = 1.0/mesh.deltaCoeffs().internalField();
    scalarField wf = mesh.weights().internalField();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const vectorField& Cf = mesh.Cf().internalField();
    const vectorField& C = mesh.C().internalField();

    #include "updateCoefficients.H"

    const tensorField& CQDWI = CQDW.internalField();
    const tensorField& CQDThetaI = CQDTheta.internalField();
    const tensorField& CMDThetaI = CMDTheta.internalField();
    const tensorField& CMDTheta2I = CMDTheta2.internalField();

    const vectorField& dQI = dQ.internalField();
    const vectorField& QI = Q.internalField();
    const vectorField& oldQI = Q.oldTime().internalField();

    const vectorField& oldMI = M.oldTime().internalField();
    const vectorField& dMI = dM.internalField();

    const vectorField& DWI = DW.internalField();
    surfaceVectorField DWf = fvc::interpolate(DW);
    const vectorField& DWfI = DWf.internalField();

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

        
        // Non-linear correction
        source[own[faceI]](0) -= dQI[faceI].x();
        source[own[faceI]](1) -= dQI[faceI].y();
        source[own[faceI]](2) -= dQI[faceI].z();

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

        // Non-linear correction
        source[nei[faceI]](0) -= -dQI[faceI].x();
        source[nei[faceI]](1) -= -dQI[faceI].y();
        source[nei[faceI]](2) -= -dQI[faceI].z();
        
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
        tensor DRP = spinTensor(Cf[faceI]-C[own[faceI]]);
        tensor DRN = spinTensor(Cf[faceI]-C[nei[faceI]]);
        tensor DRPCQ = (DRP & CQDWI[faceI]);
        tensor DRNCQ = (DRN & CQDWI[faceI]);

        // W part
        u[faceI](3,0) = DRPCQ.xx()/deltaf[faceI];
        u[faceI](3,1) = DRPCQ.xy()/deltaf[faceI];
        u[faceI](3,2) = DRPCQ.xz()/deltaf[faceI];
        
        u[faceI](4,0) = DRPCQ.yx()/deltaf[faceI];
        u[faceI](4,1) = DRPCQ.yy()/deltaf[faceI];
        u[faceI](4,2) = DRPCQ.yz()/deltaf[faceI];
        
        u[faceI](5,0) = DRPCQ.zx()/deltaf[faceI];
        u[faceI](5,1) = DRPCQ.zy()/deltaf[faceI];
        u[faceI](5,2) = DRPCQ.zz()/deltaf[faceI];

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

        // Nonlinear correction
        vector MP = (DRP & dQI[faceI]);
        source[own[faceI]](3) -= MP.x();
        source[own[faceI]](4) -= MP.y();
        source[own[faceI]](5) -= MP.z();
        
        vector MN = (DRN & dQI[faceI]);
        source[nei[faceI]](3) -= -MN.x();
        source[nei[faceI]](4) -= -MN.y();
        source[nei[faceI]](5) -= -MN.z();

        // Nonlinear correction (dW)
        vector dWxdQP = ((DWfI[faceI] - DWI[own[faceI]]) ^ QI[faceI]);
        source[own[faceI]](3) -= dWxdQP.x();
        source[own[faceI]](4) -= dWxdQP.y();
        source[own[faceI]](5) -= dWxdQP.z();

        vector dWxdQN = ((DWfI[faceI] - DWI[nei[faceI]]) ^ QI[faceI]);
        source[nei[faceI]](3) -= -dWxdQN.x();
        source[nei[faceI]](4) -= -dWxdQN.y();
        source[nei[faceI]](5) -= -dWxdQN.z();

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
        vector oldM = oldMI[faceI] + dMI[faceI];
        // vector oldM = (DLambdaI[faceI] & oldMI[faceI]);
        source[own[faceI]](3) -= oldM.x();
        source[own[faceI]](4) -= oldM.y();
        source[own[faceI]](5) -= oldM.z();
        
        source[nei[faceI]](3) -= -oldM.x();
        source[nei[faceI]](4) -= -oldM.y();
        source[nei[faceI]](5) -= -oldM.z();
    }

    if (true)
    {
        forAll (DW.boundaryField(), patchI)
        {
            const fvPatch& patch = mesh.boundary()[patchI];
            const labelList& fc = patch.faceCells();
            const vectorField& pQ = Q.boundaryField()[patchI];
            const vectorField& pDWf = DWf.boundaryField()[patchI];

            // const labelListList& cCells = mesh.cellCells();
            
            if
            (
                isA<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    DW.boundaryField()[patchI]
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
    DW.boundaryField().updateCoeffs();
    DTheta.boundaryField().updateCoeffs();
    forAll (DW.boundaryField(), patchI)
    {
        const tensorField& pCQDW = CQDW.boundaryField()[patchI];
        const tensorField& pCQDTheta = CQDTheta.boundaryField()[patchI];
        const tensorField& pCMDTheta = CMDTheta.boundaryField()[patchI];
        const tensorField& pCMDTheta2 = CMDTheta2.boundaryField()[patchI];

        const vectorField& pCf =
            mesh.Cf().boundaryField()[patchI];

        const vectorField& pDQ = dQ.boundaryField()[patchI];
        const vectorField& pOldQ = Q.oldTime().boundaryField()[patchI];

        const vectorField& pDM = dM.boundaryField()[patchI];
        const vectorField& pOldM = M.oldTime().boundaryField()[patchI];

        const fvPatch& patch = mesh.boundary()[patchI];
        const labelList& fc = patch.faceCells();

        ////// W equation

        if
        (
            isA<forceBeamDisplacementFvPatchVectorField>
            (
                DW.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementFvPatchVectorField& pDW =
                refCast<forceBeamDisplacementFvPatchVectorField>
                (
                    DW.boundaryField()[patchI]
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
                DW.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementFvPatchVectorField& pDW =
                refCast<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    DW.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh.deltaCoeffs().boundaryField()[patchI];

            vectorField tang =
            (
                DLambda.boundaryField()[patchI]
              & dRuDs.boundaryField()[patchI]
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
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
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
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMDTheta/pDelta + pCMDTheta2);//ZT

                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQDThetat[faceI]
                        & (invCM[faceI] & (pCMDTheta[faceI]/pDelta[faceI]))//ZT
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
                                  - pOldM[faceI] - pDM[faceI]
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
            
            // Add non-linear correcton of force
            vectorField pDQt = pDQ - ((tang*tang) & pDQ);
            forAll (pDQ, faceI)
            {
                source[fc[faceI]](0) -= pDQt[faceI].x();
                source[fc[faceI]](1) -= pDQt[faceI].y();
                source[fc[faceI]](2) -= pDQt[faceI].z();
            }

            // Add non-linear correcton of force (old time)
            vectorField pOldQt = pOldQ - ((tang*tang) & pOldQ);
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
                DW.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pDW =
                refCast<fixedValueFvPatchVectorField>
                (
                    DW.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh.deltaCoeffs().boundaryField()[patchI];

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
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
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
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMDTheta/pDelta + pCMDTheta2); //ZT
                
                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor CqCt =
                    (
                        pCQDTheta[faceI]
                      & (invCM[faceI] & (pCMDTheta[faceI]/pDelta[faceI]))//ZT
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
                                  - pOldM[faceI] - pDM[faceI]
                                )
                            )
                        );

                    source[fc[faceI]](0) -= thetaContrib.x();
                    source[fc[faceI]](1) -= thetaContrib.y();
                    source[fc[faceI]](2) -= thetaContrib.z();
                }
            }

            // Add non-linear correcton of force
            forAll (pDQ, faceI)
            {
                source[fc[faceI]](0) -= pDQ[faceI].x();
                source[fc[faceI]](1) -= pDQ[faceI].y();
                source[fc[faceI]](2) -= pDQ[faceI].z();
            }

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
                DTheta.boundaryField()[patchI]
            )
        )
        {
            const momentBeamRotationFvPatchVectorField& pDTheta =
                refCast<momentBeamRotationFvPatchVectorField>
                (
                    DTheta.boundaryField()[patchI]
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
                DTheta.boundaryField()[patchI]
            )
        )
        {
            const fixedValueFvPatchVectorField& pDTheta =
                refCast<fixedValueFvPatchVectorField>
                (
                    DTheta.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh.deltaCoeffs().boundaryField()[patchI];

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
                vector oldM = pOldM[faceI] + pDM[faceI];
                // vector oldM = (pDLambda[faceI] & pOldM[faceI]);
                source[fc[faceI]](3) -= oldM.x();
                source[fc[faceI]](4) -= oldM.y();
                source[fc[faceI]](5) -= oldM.z();
            }
        }

        if
        (
            isA<forceBeamDisplacementFvPatchVectorField>
            (
                DW.boundaryField()[patchI]
            )
        )
        {
            const forceBeamDisplacementFvPatchVectorField& pDW =
                refCast<forceBeamDisplacementFvPatchVectorField>
                (
                    DW.boundaryField()[patchI]
                );

            // Source contribution
            forAll (pDW, faceI)
            {
                tensor DRb =
                    spinTensor
                    (
                        pCf[faceI]
                      - C[fc[faceI]]
                    );
                vector MFb = (DRb & pDW.force()[faceI]);
              
                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }
        }
        else if
        (
            isA<axialForceTransverseDisplacementFvPatchVectorField>
            (
                DW.boundaryField()[patchI]
            )
        )
        {
            const axialForceTransverseDisplacementFvPatchVectorField& pDW =
                refCast<axialForceTransverseDisplacementFvPatchVectorField>
                (
                    DW.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh.deltaCoeffs().boundaryField()[patchI];

            vectorField tang =
            (
                DLambda.boundaryField()[patchI]
              & dRuDs.boundaryField()[patchI]
            );

            tensorField pCQDWt = pCQDW - ((tang*tang) & pCQDW);
            tensorField pCQDThetat = pCQDTheta - ((tang*tang) & pCQDTheta);

            // Diag contribution
            forAll (pDW, faceI)
            {
                tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                tensor DRbCQ = (DRb & pCQDWt[faceI]);

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
                tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                tensor DRbCQ = (DRb & pCQDWt[faceI]);

                vector curSource = (DRbCQ & pDW.refDisp()[faceI])/pDelta[faceI];

                source[fc[faceI]](3) -= curSource.x();
                source[fc[faceI]](4) -= curSource.y();
                source[fc[faceI]](5) -= curSource.z();
            }

            if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
                    );

                // vectorField pDThetat = pDTheta;
                // pDThetat -= tang*(tang & pDTheta);

                // Source contribution
                forAll (pDTheta, faceI)
                {
                    tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                    tensor DRbCQCT = (DRb & pCQDThetat[faceI]);

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
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMDTheta/pDelta+pCMDTheta2);
                                        
                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );
                    tensor CqCt =
                    (
                        pCQDThetat[faceI]
                      & (invCM[faceI] & (pCMDTheta[faceI]/pDelta[faceI]))
                    );
                    // CqCt -= (CqCt & (tang[faceI]*tang[faceI]));
                    tensor DRbCQCT = (DRb & CqCt);
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
                    tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );

                    // tensor CqCt = pCQDTheta[faceI];
                    vector thetaContrib =
                    (
                        // CqCt
                        pCQDThetat[faceI]
                      & (
                            invCM[faceI]
                          & (
                                pDTheta.moment()[faceI]
                              - pOldM[faceI] - pDM[faceI]
                            )
                        )
                    );
                    // thetaContrib -=
                    //     (thetaContrib & (tang[faceI]*tang[faceI]));
                    
                    // tensor DRbCQCT = (DRb & pCQDTheta[faceI]);
                    
                    vector MFb = (DRb & thetaContrib);
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
            forAll (pDQ, faceI)
            {
                tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                vector pQa = tang[faceI]*pDW.axialForce()[faceI];
                vector MFb = (DRb & pQa);

                source[fc[faceI]](3) -= MFb.x();
                source[fc[faceI]](4) -= MFb.y();
                source[fc[faceI]](5) -= MFb.z();
            }
            
            // Add non-linear correcton of force
            forAll (pDQ, faceI)
            {
                tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                
                vector pQt = pDQ[faceI] + pOldQ[faceI];
                pQt -= ((tang[faceI]*tang[faceI]) & pQt);
                
                vector MFb = (DRb & pQt);
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
                DW.boundaryField()[patchI]
            )
        )
        {
            Info << "fixedValueFvPatchVectorField - 2" << endl;
            
            const fixedValueFvPatchVectorField& pDW =
                refCast<fixedValueFvPatchVectorField>
                (
                    DW.boundaryField()[patchI]
                );

            scalarField pDelta =
                1.0/mesh.deltaCoeffs().boundaryField()[patchI];

            // Diag contribution
            forAll (pDW, faceI)
            {
                tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                tensor DRbCQ = (DRb & pCQDW[faceI]);

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
                tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                tensor DRbCQ = (DRb & pCQDW[faceI]);

                vector curSource = (DRbCQ & pDW[faceI])/pDelta[faceI];

                source[fc[faceI]](3) -= curSource.x();
                source[fc[faceI]](4) -= curSource.y();
                source[fc[faceI]](5) -= curSource.z();
            }

            if
            (
                isA<fixedValueFvPatchVectorField>
                (
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const fixedValueFvPatchVectorField& pDTheta =
                    refCast<fixedValueFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
                    );

                // Source contribution
                forAll (pDTheta, faceI)
                {
                    tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                    tensor DRbCQCT = (DRb & pCQDTheta[faceI]);

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
                    DTheta.boundaryField()[patchI]
                )
            )
            {
                const momentBeamRotationFvPatchVectorField& pDTheta =
                    refCast<momentBeamRotationFvPatchVectorField>
                    (
                        DTheta.boundaryField()[patchI]
                    );

                tensorField invCM = inv(pCMDTheta/pDelta+pCMDTheta2);
                    
                // Diag contribution
                forAll (pDTheta, faceI)
                {
                    tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );
                    tensor DRbCQCT =
                    (
                        DRb
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
                    tensor DRb = spinTensor( pCf[faceI] - C[fc[faceI]] );
                    tensor DRbCQCT = (DRb & pCQDTheta[faceI]);

                    vector MFb =
                        (
                            DRbCQCT
                          & (
                                invCM[faceI]
                              & (
                                    pDTheta.moment()[faceI]
                                  - pOldM[faceI]
                                  - pDM[faceI]
                                )
                            )
                        );

                    source[fc[faceI]](3) -= MFb.x();
                    source[fc[faceI]](4) -= MFb.y();
                    source[fc[faceI]](5) -= MFb.z();
                }
            }

            // Add non-linear correcton of force
            forAll (pDQ, faceI)
            {
                tensor DRb = spinTensor(pCf[faceI]-C[fc[faceI]]);
                vector MFb =
                    (DRb & (pDQ[faceI] + pOldQ[faceI]));

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

    // scalar alpha = 0.9;
    // source += (d & DWDTheta)*(1-alpha)/alpha;
    // d = d/alpha;

    
    // Block coupled solver call
    vector6 eqnRes = DWDThetaEqn.solve().initialResidual();

    maxResidual = mag(eqnRes);

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

    Info << "Max residual: " << maxResidual << endl;

    DWDThetaEqn.retrieveSolution(3, DTheta.internalField());
    DTheta.correctBoundaryConditions();

    // DTheta.relax();

    DWDThetaEqn.retrieveSolution(0, DW.internalField());
    // DW.boundaryField().updateCoeffs();
    DW.correctBoundaryConditions();

    
    // Q = (CQDW & fvc::snGrad(DW)) + (CQDW & ((II - DLambda) & dRuDs))
    //   + (DLambda & Q.oldTime());

    // DW.relax();
    // DTheta.relax();

    Q.storePrevIter();
    Q = (CQDW & fvc::snGrad(DW))
      + (CQDTheta & fvc::interpolate(DTheta)) + dQ
      + Q.oldTime();
      // + (DLambda & Q.oldTime());
    Q.relax();

    M.storePrevIter();
    M = (CMDTheta & fvc::snGrad(DTheta))
      + (CMDTheta2 & fvc::interpolate(DTheta)) + dM
      + M.oldTime();
      // + (DLambda & M.oldTime());
    M.relax();

    // Info << DW.boundaryField() << endl;
    
    DW.relax();
    DTheta.relax();
}
