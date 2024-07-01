{
    scalar ThetaResidual = GREAT;
    scalar WResidual = GREAT;

    // Initialise the block system
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

    const labelList& own = mesh().owner(); // unallocLabelList => labelList (ESI)
    const labelList& nei = mesh().neighbour();

    W_.boundaryFieldRef().updateCoeffs();
    Theta_.boundaryFieldRef().updateCoeffs();

    #include "updateCoefficients_TLNR.H"

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

    // SB added - 10/11/2023 - initial accleration and velocity values of 0th iteration
    // Newmark-beta integration scheme
    if(!steadyState() && newmark_)
    {
        if(iOuterCorr() == 0)
        {
            Accl_ = - (1/(runTime().deltaT()*betaN_))*U_.oldTime()
                - (0.5/betaN_ - 1)*Accl_.oldTime();

            U_ = U_.oldTime()
                + runTime().deltaT()*((1 - gammaN_)*Accl_.oldTime() + gammaN_*Accl_);

            dotOmega_ = - (1/(runTime().deltaT()*betaN_))*Omega_.oldTime()
                - (0.5/betaN_ - 1)*dotOmega_.oldTime();

            Omega_ = Omega_.oldTime()
                + runTime().deltaT()*((1 - gammaN_)*dotOmega_.oldTime() + gammaN_*dotOmega_);
        }
    }

    // Internal faces
    forAll(u, faceI)
    {
        //-W part (Laplacian)
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

        //- Theta part: Laplacian

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


        //- Theta part
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

        //- Laplacian part

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


        //- Theta part

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

        source[own[faceI]](3,0) -= correctedOwnExplicitMQ.x();
        source[own[faceI]](4,0) -= correctedOwnExplicitMQ.y();
        source[own[faceI]](5,0) -= correctedOwnExplicitMQ.z();

        vector correctedNeiExplicitMQ = explicitMQI[faceI];

        // Similarly here, the explicit term is not in the negative format? // understood this
        source[nei[faceI]](3,0) -= correctedNeiExplicitMQ.x();
        source[nei[faceI]](4,0) -= correctedNeiExplicitMQ.y();
        source[nei[faceI]](5,0) -= correctedNeiExplicitMQ.z();
    }

    // Boundary contributions
    forAll (W_.boundaryField(), patchI)
    {
        const fvPatch& patch = mesh().boundary()[patchI];
        const labelList& fc = patch.faceCells();

        if (patch.coupled())
        {
            notImplemented("Not ported for parallel yet!");
            //#include "updateCouplingCoeffs.H"
        }
        else
        {
                #include "applyBoundaryConditions.H"
        }
    }

    // Add distributed force
    forAll(source, cellI)
    {
        source[cellI](0,0) -= q()[cellI].x()*L()[cellI];
        source[cellI](1,0) -= q()[cellI].y()*L()[cellI];
        source[cellI](2,0) -= q()[cellI].z()*L()[cellI];
    }

    forAll(source, cellI)
    {
        source[cellI](0,0) -= rho().value()*L()[cellI]*A().value()*g().component(0).value();
        source[cellI](1,0) -= rho().value()*L()[cellI]*A().value()*g().component(1).value();
        source[cellI](2,0) -= rho().value()*L()[cellI]*A().value()*g().component(2).value();
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
    // Add body force due to gravity and buoyancy if fluid is present
    // Note that for buoyant force calculation it is assumed
    // that the entire beam is submerged in the fluid
    forAll(source, cellI)
    {
        source[cellI](0,0) -=
            (
                rho().value() - rhoFluid().value()
            )*L()[cellI]*A().value()*g().component(0).value();

        source[cellI](1,0) -=
            (
                rho().value() - rhoFluid().value()
            )*L()[cellI]*A().value()*g().component(1).value();

        source[cellI](2,0) -=
            (
                rho().value() - rhoFluid().value()
            )*L()[cellI]*A().value()*g().component(2).value();
    }

    // ground contact contribution
    if (groundContactActive_)
    {
        label cellsInContact = 0 ;
        forAll(source, cellI)
        {
            const vector coord = refW_[cellI] + W_[cellI];

            if (coord.z() < groundZ_)
            {
                cellsInContact += 1;
                source[cellI](2,0) +=
                    (2.0*gStiffness_*R()*(coord.z() - groundZ_))
                    - (2.0*gDamping_*R()*max(U_[cellI].component(2), 0));
            }
        }
        Info<< "Number of cells in contact : " << cellsInContact << endl;
    }

    // Add inertial forces
    if (!steadyState())
    {
        // First order Euler scheme
        // Add inertial force
        {
            // SB modified - (10/11/2023)
            if(newmark_)
            {
                volVectorField a = Accl_;
                vectorField QRho = ARho_*L()*a;


                forAll(source, cellI)
                {
                    source[cellI](0,0) += QRho[cellI].x();
                    source[cellI](1,0) += QRho[cellI].y();
                    source[cellI](2,0) += QRho[cellI].z();
                }

                // Add diagonal contribution (N-R method)
                scalarField QRhoCoeff =
                    -L()*ARho_/(sqr(runTime().deltaT().value())*betaN_);

                forAll(d, cellI)
                {
                    d[cellI](0,0) += QRhoCoeff[cellI];
                    d[cellI](1,1) += QRhoCoeff[cellI];
                    d[cellI](2,2) += QRhoCoeff[cellI];
                }

            }
            else
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
        }


        // Add inertial torque
        {
            if(newmark_)
            {
                // Angular acceleration
                volVectorField dotOmega = dotOmega_;

                volVectorField MRho
                    (

                        L()
                        *(
                            (RM_ & (CIRho_ & dotOmega))
                            + (RM_ & (Omega_ ^ (CIRho_ & Omega_)))
                        )
                    );
                forAll(source, cellI)
                {
                    source[cellI](3,0) += MRho[cellI].x();
                    source[cellI](4,0) += MRho[cellI].y();
                    source[cellI](5,0) += MRho[cellI].z();
                }

                // Add diagonal contribution (N-R method)- SB : Nov 2023
                // Implicit contrbution without tangent space
                volTensorField MRhoCoeff
                    (
                        L()
                        *(
                            spinTensor((RM_ & (Omega_ ^ (CIRho_ & Omega_))) + (RM_ & (CIRho_ & dotOmega)))

                            + (gammaN_/(betaN_*runTime().deltaT()))
                            *(
                                (spinTensor(RM_ & (CIRho_ & Omega_)))
                                - (RM_ & (spinTensor(Omega_) & (CIRho_ &  RM_.T())))
                            )

                            - (1/(betaN_*sqr(runTime().deltaT())))
                            *(
                                RM_ & (CIRho_ & RM_.T())
                            )
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
            else
            {
                // Angular acceleration
                volVectorField dotOmega(fvc::ddt(Omega_));

                volVectorField MRho
                    (
                        L()
                        *(
                            (RM_ & (CIRho_ & dotOmega))
                            + (RM_ & (Omega_ ^ (CIRho_ & Omega_)))
                        )
                    );

                forAll(source, cellI)
                {
                    source[cellI](3,0) += MRho[cellI].x();
                    source[cellI](4,0) += MRho[cellI].y();
                    source[cellI](5,0) += MRho[cellI].z();
                }
		        //- An acess to fluid Mesh through the parent time
                const fvMesh& fluidMesh = mesh().time().db().parent().lookupObject<fvMesh>("region0");
                 
                const vectorField beamCellCenterCoord = mesh().C() + refW_ + W_ ;

                
                //- Drag forces due to Morison's Equation
                if (dragActive_ && !steadyState())
                {
                    // Create spline using current beam points and tangents data
                    HermiteSpline spline
                    (
                        currentBeamPoints(),
                        currentBeamTangents()
                    );

                    // Evaluate dRdS - tangents to beam centreline at beam CV cell-centres
                    const vectorField& dRdScell = spline.midPointDerivatives();

                    // Tangential component of velocity vector
                    vectorField Ut
                        (
                            (
                                (U_.internalField() & dRdScell)
                                *dRdScell
                            )
                        );

                    vectorField UtHat (Ut/(mag(Ut) + SMALL));

                    // Normal component of velocity vector
                    vectorField Un
                        (
                            (
                                U_.internalField()
                                - (
                                    (U_.internalField() & dRdScell)
                                    *dRdScell
                                )
                            )
                        );

                    vectorField UnHat (Un/(mag(Un) + SMALL));
                    
                    // Scalar values of drag force (normal and tangential)
                    if (couplingActive_)
                    {
                        labelList seedCellIDs(mesh().nCells(), -1);
                        const vectorField fluidU (getFluidVelocity(fluidMesh,mesh(),beamCellCenterCoord,seedCellIDs,groundZ_,searchEngine_));
                        vectorField fluidUt
                        (
                            (
                                (fluidU & dRdScell)
                                *dRdScell
                            )
                        );
                        vectorField fluidUtHat (fluidUt/(mag(fluidUt) + SMALL));
                        // Normal component of velocity vector
                        vectorField fluidUn
                            (
                                (
                                    fluidU
                                    - (
                                        (fluidU & dRdScell)
                                        *dRdScell
                                    )
                                )
                            );
                        vectorField fluidUnHat (fluidUn/(mag(fluidUn) + SMALL));

                        const scalarField Fdn (rho().value()*Cdn_*R()*L()*((Un-fluidUn) & (Un-fluidUn)));
                        const scalarField Fdt (rho().value()*Cdt_*R()*L()*((Ut-fluidUt) & (Ut-fluidUt)));
                        forAll(source, cellI)
                        {
                            source[cellI](0,0) += Fdn[cellI]*UnHat[cellI].component(0);
                            source[cellI](1,0) += Fdn[cellI]*UnHat[cellI].component(1);
                            source[cellI](2,0) += Fdn[cellI]*UnHat[cellI].component(2);

                            source[cellI](0,0) += Fdt[cellI]*UtHat[cellI].component(0);
                            source[cellI](1,0) += Fdt[cellI]*UtHat[cellI].component(1);
                            source[cellI](2,0) += Fdt[cellI]*UtHat[cellI].component(2);
                        }

                        //todo check the force from beam to fluid => through fvOptions
                    }
                    else
                    {
                        const scalarField Fdn(rho().value()*Cdn_*R()*L()*(Un & Un));
                        const scalarField Fdt(rho().value()*Cdt_*R()*L()*(Ut & Ut));

                        forAll(source, cellI)
                        {
                            source[cellI](0,0) += Fdn[cellI]*UnHat[cellI].component(0);
                            source[cellI](1,0) += Fdn[cellI]*UnHat[cellI].component(1);
                            source[cellI](2,0) += Fdn[cellI]*UnHat[cellI].component(2);

                            source[cellI](0,0) += Fdt[cellI]*UtHat[cellI].component(0);
                            source[cellI](1,0) += Fdt[cellI]*UtHat[cellI].component(1);
                            source[cellI](2,0) += Fdt[cellI]*UtHat[cellI].component(2);
                        }

                    }
                    // Explicit drag forces included in the source vector
                    
                }
                else
                {
                    WarningIn("coupledTotalLagNewtonRaphsonBeam::evolve()")
                        << "Drag forces are zero for steady state calculation"
                        << nl
                        << "Set both 'steadyState' flag to false and "
                        << "'dragActive' flag to  true to include drag force"
                        << " contributions" << nl << endl;
                }

                // Add diagonal contribution (N-R method)-ZT code
                volTensorField MRhoCoeff
                    (
                        L()
                        *(
                            spinTensor(RM_ & (CIRho_ & Omega_))/runTime().deltaT()

                            - (RM_ & (CIRho_ & RM_.T()))/sqr(runTime().deltaT())

                            - spinTensor(RM_ & (CIRho_ & Omega_.oldTime()))/runTime().deltaT() // + sign ?

                            - spinTensor(RM_ & (spinTensor(Omega_) & (CIRho_ & Omega_)))

                            - spinTensor(RM_ & (CIRho_ & Omega_))/runTime().deltaT() // - sign

                            + (RM_ & (spinTensor(Omega_) & (CIRho_ & RM_.T())))/runTime().deltaT()
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
    }


    // Calculate equilibrium equations residual
    if (debug)
    {
        // Note: we do not user a gSum here as the beam is assumed to be on one
        // core
        const scalar eqResidual = sqrt(sum(magSqr(source)));
        Info<< "L2 norm of the equlibrium equations residual: "
            << eqResidual << endl;
    }


    // Block coupled solver call

    // Create Eigen linear solver
    BlockEigenSolverOF eigenSolver(d, l, u, own, nei);


    // Create solution vector
    Field<scalarRectangularMatrix> solVec
    (
        mesh().nCells(), scalarRectangularMatrix(6, 1, 0.0)
    );

    // Solve the linear system
    // Currently this residual is not used, Check with Seevani.
    currentResidual = eigenSolver.solve(solVec, source); // peak RAM
    Info<< "Equation Residual: " << currentResidual << endl;

    //vector6 eqnRes = WThetaEqn.solve().initialResidual();
    // vector eqnRes = vector::zero;
    // currentResidual = mag(eqnRes);


    if (iOuterCorr() == 0)
    {
        initialResidual = currentResidual;
    }

    // Copy the solution from solVec into the DW and DTheta fields
    //WThetaEqn.retrieveSolution(3, DTheta_.internalField());
    //DTheta_.boundaryField().evaluateCoupled();
    vectorField& DWI = DW_;
    vectorField& DThetaI = DTheta_;
    forAll(solVec, cellI)
    {
        DWI[cellI][vector::X] = solVec[cellI](0, 0);
        DWI[cellI][vector::Y] = solVec[cellI](1, 0);
        DWI[cellI][vector::Z] = solVec[cellI](2, 0);

        DThetaI[cellI][vector::X] = solVec[cellI](3, 0);
        DThetaI[cellI][vector::Y] = solVec[cellI](4, 0);
        DThetaI[cellI][vector::Z] = solVec[cellI](5, 0);
    }

    DW_.correctBoundaryConditions();
    DTheta_.correctBoundaryConditions();

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

    // Update displacement increment (for contact calculation of pulleys)
    WIncrement_ = W_ - W_.oldTime();

    // Update mean line velocity field
    //  SB added: (10/11/2023) Update mean line velocity field
    if (newmark_)
    {
        U_ += (1/runTime().deltaT())*(gammaN_/betaN_)*DW_;

        // SB added: (19/04/2023) Update mean line acceleration field
        Accl_ += (1/(sqr(runTime().deltaT())*betaN_))*DW_;

    }
    else
    {
        U_ = fvc::ddt(W_);
    }
 

    if (objectiveInterpolation())
    {
        // Info<< "Using objective interpolation for rotation" << endl;

        // Calculate rotation matrix correction from
        // cell-centre rotation vector correction
        volTensorField DLambda(rotationMatrix(DTheta_));

        // Update cell-centre rotation matrix
        RM_ = (DLambda & RM_);

        // Calculate mean line curvature at cell-faces
        K_ = refLambda_.T() & meanLineCurvature(RM_); // this does not work
        //K_ +=
        //    (
        //          (refLambda_.T() & Lambda_.T()) &
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
        //Info<< "Rotations are not interpolated objectively \n" << endl;
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

        if (newmark_)
        {
            // Update angular velocity
            // without tangent space
            Omega_ += (gammaN_/(betaN_*runTime().deltaT()))*(RM_.T() &  DTheta_);

            // Update angular acceleration
            dotOmega_ += (1/(betaN_*sqr(runTime().deltaT())))*(RM_.T() &  DTheta_);
        }
        else // First order Euler scheme
        {
            Omega_ = axialVector(RM_.T() & fvc::ddt(RM_));

        }
    }

    // Update axial and shear strain vector
    {
        // W_.correctBoundaryConditions();

        surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

        Gamma_ = (refLambda_.T() & ((Lambda_.T() & dRdS) - dR0Ds_));
    }

    // Q_ = (CQ_ & (Gamma_ - GammaP_));
    // M_ = (CM_ & (K_ - KP_));

    // Q_ = ((Lambda_ & refLambda_) & (CQ_ & (Gamma_ - GammaP_)));
    // M_ = ((Lambda_ & refLambda_) & (CM_ & (K_ - KP_)));

    // Calculate Q, where we ignore the orientation check
    {
        explicitQ_.setOriented(true);
        surfaceVectorField implicitQ(CQW_ & fvc::snGrad(DW_));
        implicitQ.setOriented(false);
        implicitQ += (CQTheta_ & fvc::interpolate(DTheta_));
        implicitQ.setOriented(true);
        Q_ = explicitQ_ + implicitQ;
    }

    // Calculate M, where we ignore the orientation check
    {
        explicitM_.setOriented(true);
        surfaceVectorField implicitM(CMTheta_ & fvc::snGrad(DTheta_));
        implicitM.setOriented(false);
        implicitM += (CMTheta2_ & fvc::interpolate(DTheta_));
        implicitM.setOriented(true);
        M_ = explicitM_ + implicitM;
    }

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
            // gMax
            max
            (
                mag
                (
                    Theta_.primitiveFieldRef()
                  - Theta_.oldTime().primitiveFieldRef()
                )()
            );

        if (denom < 10*SMALL)
        {
            // denom = max(gMax(mag(Theta_.internalField())), SMALL);
            denom = 1.0;
        }

        ThetaResidual =
            // gMax(mag(DTheta_.internalField())());
            max(mag(DTheta_.primitiveFieldRef())());
        // ThetaResidual =
        //     gMax(mag(DTheta_.internalField()))/denom;
    }

    // Calculate DW residual
    {
        scalar denom =
            max //gMax
            (
                mag
                (
                    W_.primitiveFieldRef()
                  - W_.oldTime().primitiveFieldRef()
                )()
            );

        if (denom < 10*SMALL)
        {
            // denom = max(gMax(mag(W_.internalField())), SMALL);
            denom = 1.0;
        }

        WResidual =
            // gMax(mag(DW_.internalField())());
            max(mag(DW_.primitiveFieldRef())());
        // WResidual =
        //     gMax(mag(DW_.internalField()))/denom;
    }

    if (debug)
    {
        Info<< "Theta residual: " << ThetaResidual << endl;
        Info<< "W residual: " << WResidual << endl;
    }

    // this currentResidual is not normalized, check with Seevani
    currentResidual = max(WResidual, ThetaResidual);

    if (iOuterCorr() == 0)
    {
        initialResidual = currentResidual;
    }
}
