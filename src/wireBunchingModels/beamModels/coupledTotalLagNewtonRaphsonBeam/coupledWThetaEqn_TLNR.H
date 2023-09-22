{
    // Info << "Using coupled solver for DW and DTheta" << endl;

    scalar ThetaResidual = GREAT;
    scalar WResidual = GREAT;

    // Theta_.storePrevIter();
    // W_.storePrevIter();
    // Lambda_.storePrevIter();
    // Gamma_.storePrevIter();
    // K_.storePrevIter();

    // Initialise the block system
    // multibeamFvBlockMatrix WThetaEqn(WTheta_, *this);
    // fvBlockMatrix<vector6> WThetaEqn(WTheta_);

    // Pout << WTheta_.boundaryField() << endl;

    // Grab block diagonal and set it to zero
    //tensor6Field& d = WThetaEqn.diag().asSquare();
    //d = tensor6::zero;
    //scalarSquareMatrix myMatrix(6, 6, 0.0); // 6 x 6 matrix with 0.0 values

    Field<scalarSquareMatrix> d(mesh().nCells(), scalarSquareMatrix(6, 0.0));

    // Grab off-diagonal and set it to zero
    Field<scalarSquareMatrix> u(mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0));
    // tensor6Field& u = WThetaEqn.upper().asSquare();
    // u = tensor6::zero;

    // Grab off-diagonal and set it to zero
    Field<scalarSquareMatrix> l(mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0));
    // tensor6Field& l = WThetaEqn.lower().asSquare();
    // l = tensor6::zero;

    // Grap source and set it to zero
    Field<scalarRectangularMatrix> source
    (
        mesh().nCells(), scalarRectangularMatrix(6, 1, 0.0)
    );
    // vector6Field& source = WThetaEqn.source();
    // source = vector6::zero;

    scalarField deltaf = 1.0/mesh().deltaCoeffs().internalField();
    const scalarField& wf = mesh().weights().internalField();

    // const vectorField& dR0DsI = dR0Ds_.internalField();
    // const vectorField& iI = i_.internalField();

    const labelList& own = mesh().owner(); // unallocLabelList => labelList (ESI)
    const labelList& nei = mesh().neighbour();

    W_.boundaryFieldRef().updateCoeffs();
    Theta_.boundaryFieldRef().updateCoeffs();

    #include "updateCoefficients_TLNR.H"

    // Info << "explicitQ_ \n " << explicitQ_ << endl;

    const tensorField& CQWI = CQW_.internalField();
    const tensorField& CQDThetaI = CQDTheta_.internalField();
    const tensorField& CQThetaI = CQTheta_.internalField();
    const tensorField& CMThetaI = CMTheta_.internalField();
    const tensorField& CMTheta2I = CMTheta2_.internalField();
    const tensorField& CMQWI = CMQW_.internalField();
    const tensorField& CMQThetaI = CMQTheta_.internalField();

    const vectorField& explicitQI = explicitQ_.internalField();
    const vectorField& explicitMI = explicitM_.internalField();
    const vectorField& explicitMQI = explicitMQ_.internalField();

    // const vectorField& dRdSI = dRdS_.internalField();
    // const tensorField& RMI = RM_.internalField();
    // const tensorField& refRMI = refRM_.internalField();

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

        // Theta part: Laplacian

        u[faceI](0,3) += CQDThetaI[faceI].xx()/deltaf[faceI];
        u[faceI](0,4) += CQDThetaI[faceI].xy()/deltaf[faceI];
        u[faceI](0,5) += CQDThetaI[faceI].xz()/deltaf[faceI];

        u[faceI](1,3) += CQDThetaI[faceI].yx()/deltaf[faceI];
        u[faceI](1,4) += CQDThetaI[faceI].yy()/deltaf[faceI];
        u[faceI](1,5) += CQDThetaI[faceI].yz()/deltaf[faceI];

        u[faceI](2,3) += CQDThetaI[faceI].zx()/deltaf[faceI];
        u[faceI](2,4) += CQDThetaI[faceI].zy()/deltaf[faceI];
        u[faceI](2,5) += CQDThetaI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](0,3) += -CQDThetaI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](0,4) += -CQDThetaI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](0,5) += -CQDThetaI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](1,3) += -CQDThetaI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](1,4) += -CQDThetaI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](1,5) += -CQDThetaI[faceI].yz()/deltaf[faceI];

        d[own[faceI]](2,3) += -CQDThetaI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](2,4) += -CQDThetaI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](2,5) += -CQDThetaI[faceI].zz()/deltaf[faceI];

        l[faceI](0,3) += CQDThetaI[faceI].xx()/deltaf[faceI];
        l[faceI](0,4) += CQDThetaI[faceI].xy()/deltaf[faceI];
        l[faceI](0,5) += CQDThetaI[faceI].xz()/deltaf[faceI];

        l[faceI](1,3) += CQDThetaI[faceI].yx()/deltaf[faceI];
        l[faceI](1,4) += CQDThetaI[faceI].yy()/deltaf[faceI];
        l[faceI](1,5) += CQDThetaI[faceI].yz()/deltaf[faceI];

        l[faceI](2,3) += CQDThetaI[faceI].zx()/deltaf[faceI];
        l[faceI](2,4) += CQDThetaI[faceI].zy()/deltaf[faceI];
        l[faceI](2,5) += CQDThetaI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](0,3) += -CQDThetaI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](0,4) += -CQDThetaI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](0,5) += -CQDThetaI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](1,3) += -CQDThetaI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](1,4) += -CQDThetaI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](1,5) += -CQDThetaI[faceI].yz()/deltaf[faceI];

        d[nei[faceI]](2,3) += -CQDThetaI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](2,4) += -CQDThetaI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](2,5) += -CQDThetaI[faceI].zz()/deltaf[faceI];


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

        source[own[faceI]](0,0) -= explicitQI[faceI].x();
        source[own[faceI]](1,0) -= explicitQI[faceI].y();
        source[own[faceI]](2,0) -= explicitQI[faceI].z();

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

        source[nei[faceI]](0,0) -= -explicitQI[faceI].x();
        source[nei[faceI]](1,0) -= -explicitQI[faceI].y();
        source[nei[faceI]](2,0) -= -explicitQI[faceI].z();


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


        // Theta part

        u[faceI](3,3) += (1-wf[faceI])*CMTheta2I[faceI].xx();
        u[faceI](3,4) += (1-wf[faceI])*CMTheta2I[faceI].xy();
        u[faceI](3,5) += (1-wf[faceI])*CMTheta2I[faceI].xz();

        u[faceI](4,3) += (1-wf[faceI])*CMTheta2I[faceI].yx();
        u[faceI](4,4) += (1-wf[faceI])*CMTheta2I[faceI].yy();
        u[faceI](4,5) += (1-wf[faceI])*CMTheta2I[faceI].yz();

        u[faceI](5,3) += (1-wf[faceI])*CMTheta2I[faceI].zx();
        u[faceI](5,4) += (1-wf[faceI])*CMTheta2I[faceI].zy();
        u[faceI](5,5) += (1-wf[faceI])*CMTheta2I[faceI].zz();

        d[own[faceI]](3,3) += wf[faceI]*CMTheta2I[faceI].xx();
        d[own[faceI]](3,4) += wf[faceI]*CMTheta2I[faceI].xy();
        d[own[faceI]](3,5) += wf[faceI]*CMTheta2I[faceI].xz();

        d[own[faceI]](4,3) += wf[faceI]*CMTheta2I[faceI].yx();
        d[own[faceI]](4,4) += wf[faceI]*CMTheta2I[faceI].yy();
        d[own[faceI]](4,5) += wf[faceI]*CMTheta2I[faceI].yz();

        d[own[faceI]](5,3) += wf[faceI]*CMTheta2I[faceI].zx();
        d[own[faceI]](5,4) += wf[faceI]*CMTheta2I[faceI].zy();
        d[own[faceI]](5,5) += wf[faceI]*CMTheta2I[faceI].zz();


        l[faceI](3,3) += -wf[faceI]*CMTheta2I[faceI].xx();
        l[faceI](3,4) += -wf[faceI]*CMTheta2I[faceI].xy();
        l[faceI](3,5) += -wf[faceI]*CMTheta2I[faceI].xz();

        l[faceI](4,3) += -wf[faceI]*CMTheta2I[faceI].yx();
        l[faceI](4,4) += -wf[faceI]*CMTheta2I[faceI].yy();
        l[faceI](4,5) += -wf[faceI]*CMTheta2I[faceI].yz();

        l[faceI](5,3) += -wf[faceI]*CMTheta2I[faceI].zx();
        l[faceI](5,4) += -wf[faceI]*CMTheta2I[faceI].zy();
        l[faceI](5,5) += -wf[faceI]*CMTheta2I[faceI].zz();

        d[nei[faceI]](3,3) += -(1-wf[faceI])*CMTheta2I[faceI].xx();
        d[nei[faceI]](3,4) += -(1-wf[faceI])*CMTheta2I[faceI].xy();
        d[nei[faceI]](3,5) += -(1-wf[faceI])*CMTheta2I[faceI].xz();

        d[nei[faceI]](4,3) += -(1-wf[faceI])*CMTheta2I[faceI].yx();
        d[nei[faceI]](4,4) += -(1-wf[faceI])*CMTheta2I[faceI].yy();
        d[nei[faceI]](4,5) += -(1-wf[faceI])*CMTheta2I[faceI].yz();

        d[nei[faceI]](5,3) += -(1-wf[faceI])*CMTheta2I[faceI].zx();
        d[nei[faceI]](5,4) += -(1-wf[faceI])*CMTheta2I[faceI].zy();
        d[nei[faceI]](5,5) += -(1-wf[faceI])*CMTheta2I[faceI].zz();

        // Explicit part
        source[own[faceI]](3,0) -= explicitMI[faceI].x();
        source[own[faceI]](4,0) -= explicitMI[faceI].y();
        source[own[faceI]](5,0) -= explicitMI[faceI].z();

        source[nei[faceI]](3,0) -= -explicitMI[faceI].x();
        source[nei[faceI]](4,0) -= -explicitMI[faceI].y();
        source[nei[faceI]](5,0) -= -explicitMI[faceI].z();



        //---- (dr x  Q) term

        // W part
        u[faceI](3,0) += CMQWI[faceI].xx()/deltaf[faceI];
        u[faceI](3,1) += CMQWI[faceI].xy()/deltaf[faceI];
        u[faceI](3,2) += CMQWI[faceI].xz()/deltaf[faceI];

        u[faceI](4,0) += CMQWI[faceI].yx()/deltaf[faceI];
        u[faceI](4,1) += CMQWI[faceI].yy()/deltaf[faceI];
        u[faceI](4,2) += CMQWI[faceI].yz()/deltaf[faceI];

        u[faceI](5,0) += CMQWI[faceI].zx()/deltaf[faceI];
        u[faceI](5,1) += CMQWI[faceI].zy()/deltaf[faceI];
        u[faceI](5,2) += CMQWI[faceI].zz()/deltaf[faceI];

        d[own[faceI]](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        d[own[faceI]](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        d[own[faceI]](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        d[own[faceI]](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        d[own[faceI]](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        d[own[faceI]](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        d[own[faceI]](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        d[own[faceI]](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        d[own[faceI]](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];

        l[faceI](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        l[faceI](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        l[faceI](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        l[faceI](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        l[faceI](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        l[faceI](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        l[faceI](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        l[faceI](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        l[faceI](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];

        d[nei[faceI]](3,0) += CMQWI[faceI].xx()/deltaf[faceI];
        d[nei[faceI]](3,1) += CMQWI[faceI].xy()/deltaf[faceI];
        d[nei[faceI]](3,2) += CMQWI[faceI].xz()/deltaf[faceI];

        d[nei[faceI]](4,0) += CMQWI[faceI].yx()/deltaf[faceI];
        d[nei[faceI]](4,1) += CMQWI[faceI].yy()/deltaf[faceI];
        d[nei[faceI]](4,2) += CMQWI[faceI].yz()/deltaf[faceI];

        d[nei[faceI]](5,0) += CMQWI[faceI].zx()/deltaf[faceI];
        d[nei[faceI]](5,1) += CMQWI[faceI].zy()/deltaf[faceI];
        d[nei[faceI]](5,2) += CMQWI[faceI].zz()/deltaf[faceI];

        // Theta part
        u[faceI](3,3) += (1-wf[faceI])*CMQThetaI[faceI].xx();
        u[faceI](3,4) += (1-wf[faceI])*CMQThetaI[faceI].xy();
        u[faceI](3,5) += (1-wf[faceI])*CMQThetaI[faceI].xz();

        u[faceI](4,3) += (1-wf[faceI])*CMQThetaI[faceI].yx();
        u[faceI](4,4) += (1-wf[faceI])*CMQThetaI[faceI].yy();
        u[faceI](4,5) += (1-wf[faceI])*CMQThetaI[faceI].yz();

        u[faceI](5,3) += (1-wf[faceI])*CMQThetaI[faceI].zx();
        u[faceI](5,4) += (1-wf[faceI])*CMQThetaI[faceI].zy();
        u[faceI](5,5) += (1-wf[faceI])*CMQThetaI[faceI].zz();

        d[own[faceI]](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        d[own[faceI]](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        d[own[faceI]](3,5) += wf[faceI]*CMQThetaI[faceI].xz();

        d[own[faceI]](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        d[own[faceI]](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        d[own[faceI]](4,5) += wf[faceI]*CMQThetaI[faceI].yz();

        d[own[faceI]](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        d[own[faceI]](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        d[own[faceI]](5,5) += wf[faceI]*CMQThetaI[faceI].zz();


		// Shouldn't the lower and diag[nei] coefficients have a negative sign?
        l[faceI](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        l[faceI](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        l[faceI](3,5) += wf[faceI]*CMQThetaI[faceI].xz();

        l[faceI](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        l[faceI](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        l[faceI](4,5) += wf[faceI]*CMQThetaI[faceI].yz();

        l[faceI](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        l[faceI](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        l[faceI](5,5) += wf[faceI]*CMQThetaI[faceI].zz();

        d[nei[faceI]](3,3) += (1-wf[faceI])*CMQThetaI[faceI].xx();
        d[nei[faceI]](3,4) += (1-wf[faceI])*CMQThetaI[faceI].xy();
        d[nei[faceI]](3,5) += (1-wf[faceI])*CMQThetaI[faceI].xz();

        d[nei[faceI]](4,3) += (1-wf[faceI])*CMQThetaI[faceI].yx();
        d[nei[faceI]](4,4) += (1-wf[faceI])*CMQThetaI[faceI].yy();
        d[nei[faceI]](4,5) += (1-wf[faceI])*CMQThetaI[faceI].yz();

        d[nei[faceI]](5,3) += (1-wf[faceI])*CMQThetaI[faceI].zx();
        d[nei[faceI]](5,4) += (1-wf[faceI])*CMQThetaI[faceI].zy();
        d[nei[faceI]](5,5) += (1-wf[faceI])*CMQThetaI[faceI].zz();

        // Explicit part
        vector correctedOwnExplicitMQ = explicitMQI[faceI];

        // if (false)
        // {
        //     correctedOwnExplicitMQ =
        //         (
        //             refRMI[own[faceI]].T()
        //           & (RMI[own[faceI]].T() & correctedOwnExplicitMQ)
        //         );

        //     correctedOwnExplicitMQ.x() = 0;

        //     correctedOwnExplicitMQ =
        //         (
        //             RMI[own[faceI]]
        //           & (refRMI[own[faceI]] & correctedOwnExplicitMQ)
        //         );
        // }

        source[own[faceI]](3,0) -= correctedOwnExplicitMQ.x();
        source[own[faceI]](4,0) -= correctedOwnExplicitMQ.y();
        source[own[faceI]](5,0) -= correctedOwnExplicitMQ.z();

        vector correctedNeiExplicitMQ = explicitMQI[faceI];

        // correctedNeiExplicitMQ =
        //     (
        //         refRMI[nei[faceI]].T()
        //       & (RMI[nei[faceI]].T() & correctedNeiExplicitMQ)
        //     );
        // correctedNeiExplicitMQ.x() = 0;
        // correctedNeiExplicitMQ =
        //     (RMI[nei[faceI]] & (refRMI[nei[faceI]] & correctedNeiExplicitMQ));

		// Similarly here, the explicit term is not in the negative format? // understood this
        source[nei[faceI]](3,0) -= correctedNeiExplicitMQ.x();
        source[nei[faceI]](4,0) -= correctedNeiExplicitMQ.y();
        source[nei[faceI]](5,0) -= correctedNeiExplicitMQ.z();
    }

    // Boundary contributions
    forAll (W_.boundaryField(), patchI)
    {
        const fvPatch& patch = mesh().boundary()[patchI];
        //const labelList& fc = patch.faceCells();

        if (patch.coupled())
        {
            notImplemented("Not ported for parallel yet!");
            //#include "updateCouplingCoeffs.H"
        }
        else
        {
            // Philip: we do need this, but you can comment it initially
            notImplemented("We need to port applyBoundaryConditions.H!");
            //#include "applyBoundaryConditions.H"
        }
    }

    // Add distributed force
    forAll(source, cellI)
    {
        source[cellI](0,0) -= q()[cellI].x()*L()[cellI];
        source[cellI](1,0) -= q()[cellI].y()*L()[cellI];
        source[cellI](2,0) -= q()[cellI].z()*L()[cellI];
    }

    // Add point forces
    forAll(pointForces(), pfI)
    {
        // Get beam relative coordinates
        label cellI = pointForces()[pfI].first().first();
        scalar zeta = pointForces()[pfI].first().second();

        vector F0 = pointForces()[pfI].second()(runTime().value());

        source[cellI](0,0) -= F0.x();
        source[cellI](1,0) -= F0.y();
        source[cellI](2,0) -= F0.z();

        surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

        const surfaceScalarField& dc = mesh().deltaCoeffs();

        vector DR = vector::zero;

        if (zeta > SMALL)
        {
            label faceID = own.find(cellI);
            if (faceID == -1) // last cell
            {
                const labelList& faceCells =
                    mesh().boundary()[endPatchIndex()].faceCells();

                label bFaceID = faceCells.find(cellI);

                DR = zeta*dRdS.boundaryField()[endPatchIndex()][bFaceID]
                   /dc.boundaryField()[endPatchIndex()][bFaceID];
            }
            else
            {
                DR = 0.5*zeta*dRdS.internalField()[faceID]*deltaf[faceID];
            }
        }
        else if (zeta < -SMALL)
        {
            label faceID = nei.find(cellI);
            if (faceID == -1) // first cell
            {
                const labelList& faceCells =
                    mesh().boundary()[startPatchIndex()].faceCells();

                label bFaceID = faceCells.find(cellI);

                DR = zeta*dRdS.boundaryField()[startPatchIndex()][bFaceID]
                   /dc.boundaryField()[startPatchIndex()][bFaceID];
            }
            else
            {
                DR = 0.5*zeta*dRdS.internalField()[faceID]*deltaf[faceID];
            }
        }

        vector M0 = (spinTensor(DR) & F0);

        source[cellI](3,0) -= M0.x();
        source[cellI](4,0) -= M0.y();
        source[cellI](5,0) -= M0.z();
    }

    // Add distributed moment
    forAll(source, cellI)
    {
        source[cellI](3,0) -= m()[cellI].x()*L()[cellI];
        source[cellI](4,0) -= m()[cellI].y()*L()[cellI];
        source[cellI](5,0) -= m()[cellI].z()*L()[cellI];
    }

    // Add inertial forces
    if (!steadyState())
    {
        // First order Euler scheme

        // Add inertial force
        {
            volVectorField a(fvc::ddt(U_));
            vectorField QRho = rho().value()*A().value()*L()*a;

            forAll(source, cellI)
            {
                source[cellI](0,0) += QRho[cellI].x();
                source[cellI](1,0) += QRho[cellI].y();
                source[cellI](2,0) += QRho[cellI].z();
            }

            // Add diagonal contribution (N-R method)
            scalarField QRhoCoeff =
                -L()*rho().value()*A().value()/sqr(runTime().deltaT().value());

            forAll(d, cellI)
            {
                d[cellI](0,0) += QRhoCoeff[cellI];
                d[cellI](1,1) += QRhoCoeff[cellI];
                d[cellI](2,2) += QRhoCoeff[cellI];
            }
        }


        // Add inertial torque
        {
            volVectorField dotOmega(fvc::ddt(Omega_));
            volVectorField MRho
            (
                L()*rho()
               *(
                    (RM_ & (CI_ & dotOmega))
                  + (RM_ & (Omega_ ^ (CI_ & Omega_)))
                )
            );

            forAll(source, cellI)
            {
                source[cellI](3,0) += MRho[cellI].x();
                source[cellI](4,0) += MRho[cellI].y();
                source[cellI](5,0) += MRho[cellI].z();
            }

            // Add diagonal contribution (N-R method)
            volTensorField MRhoCoeff
            (
                L()*rho()
               *(
                    spinTensor(RM_ & (CI_ & Omega_))/runTime().deltaT()

                  - (RM_ & (CI_ & RM_.T()))/sqr(runTime().deltaT())

                    - spinTensor(RM_ & (CI_ & Omega_.oldTime()))/runTime().deltaT() // + sign ?

                  - spinTensor(RM_ & (spinTensor(Omega_) & (CI_ & Omega_)))

                    - spinTensor(RM_ & (CI_ & Omega_))/runTime().deltaT() // - sign

                  + (RM_ & (spinTensor(Omega_) & (CI_ & RM_.T())))/runTime().deltaT()
               )
            );

            forAll(d, cellI)
            {
                d[cellI](3,3) += MRhoCoeff[cellI].xx();
                d[cellI](3,4) += MRhoCoeff[cellI].xy();
                d[cellI](3,5) += MRhoCoeff[cellI].xz();

                d[cellI](4,3) += MRhoCoeff[cellI].yx();
                d[cellI](4,4) += MRhoCoeff[cellI].yy();
                d[cellI](4,5) += MRhoCoeff[cellI].yz();

                d[cellI](5,3) += MRhoCoeff[cellI].zx();
                d[cellI](5,4) += MRhoCoeff[cellI].zy();
                d[cellI](5,5) += MRhoCoeff[cellI].zz();
            }
        }
    }

    // Apply beam-to-beam line contact
    //applyPointContact(WThetaEqn);

    // Apply beam-to-beam line contact
    //applyLineContact(WThetaEqn);

    // Apply conical pulley contact
    // applyConicalPulleysContactNew(WThetaEqn);

    // Apply conical pulley contact
    // applyToroidalPulleysContact(WThetaEqn);


    //scalar nCV = mesh().cellZones().size();
    // Calculate equilibrium equations residual
    if (debug)
    {
        scalar eqResidual = sqrt(gSum(magSqr(source))); // magSqr ?
        Info << "L2 norm of the equlibrium equations residual: "
             << eqResidual << endl;
    }


    // Block coupled solver call
    // PhilipC: we need to replace this with our own solver, e.g. Eigen or PETSc
    notImplemented("We need to implement out own linear solver, e.g. Eigen or PETSc");
    //vector6 eqnRes = WThetaEqn.solve().initialResidual();
    vector eqnRes = vector::zero;
    currentResidual = mag(eqnRes);


    if (iOuterCorr() == 0)
    {
        initialResidual = currentResidual;
    }

    //WThetaEqn.retrieveSolution(3, DTheta_.internalField());
    //DTheta_.boundaryField().evaluateCoupled();
    notImplemented("We need to implement Retrieve solution");
    // DTheta_.correctBoundaryConditions();

    // DTheta_.internalField().replace(0, 0);

    Theta_.primitiveFieldRef() +=
        (RM_.internalField().T() & DTheta_.internalField()); // Get back to this

    // Theta_.internalField().replace(0, 0);

    Theta_.correctBoundaryConditions();
    Theta_.storePrevIter();

    //WThetaEqn.retrieveSolution(0, DW_.internalField());
    //DW_.boundaryField().evaluateCoupled();
    W_.primitiveFieldRef() += DW_.internalField();

    W_.correctBoundaryConditions();
    W_.storePrevIter();

    // Info << DW_.internalField() << endl;

    // Update displacement increment (for contact calculation)
    WIncrement_ = W_ - W_.oldTime();

    // Update mean line velocity field
    U_ = fvc::ddt(W_);

    // Retrieve lagrange multipliers
    //WThetaEqn.retriveLagrangeMultipliers();


    if (objectiveInterpolation())
    {
        // Info << "Using objective interpolation for rotation" << endl;

        // Calculate rotation matrix correction from
        // cell-centre rotation vector correction
        volTensorField DLambda(rotationMatrix(DTheta_));

        // Update cell-centre rotation matrix
        RM_ = (DLambda & RM_);

        // Calculate mean line curvature at cell-faces
         K_ = refLambda_.T() & meanLineCurvature(RM_); // this does not work
        //K_ +=
        //    (
        // 	    (refLambda_.T() & Lambda_.T()) &
        //        meanLineCurvature(DLambda)
        //    );

        // Objective cell-to-face interpolation of rotation matrix correction
        //surfaceTensorField DLambdaf =
         //   interpolateRotationMatrix(DLambda);

        // Update cell-face rotation matrix
        //Lambda_ = (DLambdaf & Lambda_);
         Lambda_ = interpolateRotationMatrix(RM_); // this does not work
    }
    else
    {
    	Info << "Rotations are not interpolated objectively \n" << endl;
        surfaceVectorField DThetaf(fvc::interpolate(DTheta_));

        surfaceScalarField magDThetaf(mag(DThetaf) + SMALL);
        surfaceTensorField DThetaHat(spinTensor(DThetaf));

        dimensionedTensor I("I", dimless, tensor::I);

        // Tangent operator
        surfaceTensorField DT
        (
            (Foam::sin(magDThetaf)/magDThetaf)*I
          + ((1.0-Foam::sin(magDThetaf)/magDThetaf)/sqr(magDThetaf))
           *(DThetaf*DThetaf)
          + ((1.0-Foam::cos(magDThetaf))/sqr(magDThetaf))*DThetaHat
        );
        // Update bending strain vector
        K_ +=
            (
                (refLambda_.T() & Lambda_.T())
              // & fvc::snGrad(DTheta_)
              & (DT.T() & fvc::snGrad(DTheta_))
            );

        // Rodrigues formula
        surfaceTensorField DLambdaf(rotationMatrix(DThetaf));


        // Update rotation matrix
        Lambda_ = (DLambdaf & Lambda_);


        // Update cell-centre rotation matrix
        volTensorField DLambda(rotationMatrix(DTheta_));
        RM_ = (DLambda & RM_);
        // interpolateRotationMatrix(*this, Lambda_, RM_);


        // Update angular velocity
        Omega_ = axialVector(RM_.T() & fvc::ddt(RM_));
    }

    // Update axial and shear strain vector
    {
        // W_.correctBoundaryConditions();

        surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

        // Gamma_ =
        //     (
        //         (kirchhoffTransTensor_ & (refLambda_.T() & Lambda_.T()))
        //       & dRdS
        //     )
        //   - (refLambda_.T() & dR0Ds_);

        Gamma_ = (refLambda_.T() & ((Lambda_.T() & dRdS) - dR0Ds_));
    }

    // Q_ = (CQ_ & (Gamma_ - GammaP_));
    // M_ = (CM_ & (K_ - KP_));

    // Q_ = ((Lambda_ & refLambda_) & (CQ_ & (Gamma_ - GammaP_)));
    // M_ = ((Lambda_ & refLambda_) & (CM_ & (K_ - KP_)));

    // if (updatedLagrangian_)
    // {
    //    // Q_ += (Lambda_ & Q_.oldTime());
    //    M_ += (Lambda_ & M_.oldTime());
    // }


  //  Info << "explicitQ_ \n " << explicitQ_ << endl;

    Q_ = explicitQ_ +
        (CQW_ & fvc::snGrad(DW_)) +
        (CQTheta_ & fvc::interpolate(DTheta_));

  //   Info << "Q: \n" << Q_ << endl;

    M_ = explicitM_ +
        (CMTheta_ & fvc::snGrad(DTheta_)) +
        (CMTheta2_ & fvc::interpolate(DTheta_));

    // Info << "M: \n" << M << endl;

    // surfaceVectorField rPrimeXf = //explicitMQ_;
    //     explicitMQ_ +
    //     (CMQW_ & fvc::snGrad(DW_)) +
    //     (CMQTheta_ & fvc::interpolate(DTheta_));

   //  Info << "rPrimeXf: \n" << rPrimeXf << endl;

    // Calculate axial force
    {
        // surfaceVectorField t = (Lambda_ & dR0Ds_);

        // dRdS /= mag(dRdS);
        // Qa_ = (dRdS & Q_);

        Qa_ = (Lambda_.T() & Q_)().component(0);
    }

    // Calculate DTheta residual
    {
        scalar denom =
            gMax
            (
                mag
                (
                    Theta_.internalField()
                  - Theta_.oldTime().internalField()
                )()
            );

        if (denom < 10*SMALL)
        {
            // denom = max(gMax(mag(Theta_.internalField())), SMALL);
            denom = 1.0;
        }

        ThetaResidual =
            gMax(mag(DTheta_.internalField())());
        // ThetaResidual =
        //     gMax(mag(DTheta_.internalField()))/denom;
    }

    // Calculate DW residual
    {
        scalar denom =
            gMax
            (
                mag
                (
                    W_.internalField()
                  - W_.oldTime().internalField()
                )()
            );

        if (denom < 10*SMALL)
        {
            // denom = max(gMax(mag(W_.internalField())), SMALL);
            denom = 1.0;
        }

        WResidual =
            gMax(mag(DW_.internalField())());
        // WResidual =
        //     gMax(mag(DW_.internalField()))/denom;
    }

    if (debug)
    {
        Info << "Theta residual: " << ThetaResidual << endl;
        Info << "W residual: " << WResidual << endl;
    }


	if(false)
	{
		//Info << "Iter: " << iOuterCorr() << endl;
		Info << "Displacement" << W_.boundaryField() << endl;
	}
    currentResidual = max(WResidual, ThetaResidual);

    if (iOuterCorr() == 0)
    {
        initialResidual = currentResidual;
    }

    // if (plasticityStressResultantReturnPtr_.valid())
    // {
    //     GammaP_.storePrevIter();
    //     KP_.storePrevIter();

    //     if (!completedElasticPrediction)
    //     {
    //         currentMaterialResidual = 0;

    //         if
    //         (
    //             iOuterCorr()
    //           > plasticityStressResultantReturnPtr_().nElPredCorr()
    //          ||
    //             currentResidual
    //           < plasticityStressResultantReturnPtr_().elPredTol()
    //         )
    //         {
    //             Info << "nElPredictors: " << iOuterCorr() << endl;
    //             Info << "currResidual: " << currentResidual << endl;
    //             completedElasticPrediction = true;
    //         }
    //     }

    //     if (completedElasticPrediction)
    //     {
    //         if (debug)
    //         {
    //             Info << "Update plasticity" << endl;
    //         }
    //         plasticityStressResultantReturnPtr_().correct();

    //         GammaP_ = GammaP_ + plasticityStressResultantReturnPtr_().DGammaP();
    //         KP_ = KP_ + plasticityStressResultantReturnPtr_().DKP();

    //         // Calculate material residual
    //         {
    //             // scalar GammaPResidual =
    //             //     max
    //             //     (
    //             //         mag
    //             //         (
    //             //             GammaP_
    //             //           - GammaP_.prevIter()
    //             //         )()
    //             //     ).value()
    //             //    /max
    //             //     (
    //             //         (
    //             //             dimensionedScalar
    //             //             (
    //             //                 "SMALL",
    //             //                 GammaP_.dimensions(),
    //             //                 SMALL
    //             //             )
    //             //           + mag(GammaP_.prevIter())
    //             //         )()
    //             //      ).value();

    //             scalar denom =
    //                 gMax
    //                 (
    //                     mag
    //                     (
    //                         KP_
    //                       - KP_.oldTime()
    //                     )()
    //                 ); //.value();

    //             if (denom < 10*SMALL)
    //             {
    //                 // denom = max(gMax(mag(W_.internalField())), SMALL);
    //                 denom = 1.0;
    //             }

    //             scalar KPResidual =
    //                 gMax
    //                 (
    //                     mag
    //                     (
    //                         KP_
    //                       - KP_.prevIter()
    //                     )()
    //                 );//.value();
    //                 // max
    //                 // (
    //                 //     (
    //                 //         dimensionedScalar
    //                 //         (
    //                 //             "SMALL",
    //                 //             KP_.dimensions(),
    //                 //             SMALL
    //                 //         )
    //                 //       + mag(KP_.prevIter())
    //                 //     )()
    //                 // ).value();

    //             currentMaterialResidual = KPResidual;
    //             // currentMaterialResidual =
    //             //     max(GammaPResidual, KPResidual);

    //             if (debug)
    //             {
    //                 Info << "Material residual: "
    //                      << currentMaterialResidual << endl;
    //             }
    //         }

    //         CDQ_ = plasticityStressResultantReturnPtr_().DQDGamma();
    //         CDM_ = plasticityStressResultantReturnPtr_().DMDK();

    //         CDMDGamma_ = plasticityStressResultantReturnPtr_().DMDGamma();
    //         CDQDK_ = plasticityStressResultantReturnPtr_().DQDK();

    //         GammaP_.relax();
    //         KP_.relax();
    //     }
    // }
}
