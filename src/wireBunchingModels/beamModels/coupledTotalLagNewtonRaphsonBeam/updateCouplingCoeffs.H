{
    // Pout << "Update coupling coeffs for patch: " << patch.name() << endl;
    
    tensor6Field& cu = WThetaEqn.coupleUpper()[patchI].asSquare();
    cu = tensor6::zero;
    
    tensor6Field& cl = WThetaEqn.coupleLower()[patchI].asSquare();
    cl = tensor6::zero;

    scalarField deltaf = 1.0/mesh().deltaCoeffs().boundaryField()[patchI];
    const scalarField& wf = mesh().weights().boundaryField()[patchI];

    const tensorField& CQWI = CQW_.boundaryField()[patchI];
    const tensorField& CQThetaI = CQTheta_.boundaryField()[patchI];
    const tensorField& CMThetaI = CMTheta_.boundaryField()[patchI];
    const tensorField& CMTheta2I = CMTheta2_.boundaryField()[patchI];
    const tensorField& CMQWI = CMQW_.boundaryField()[patchI];
    const tensorField& CMQThetaI = CMQTheta_.boundaryField()[patchI];

    const vectorField& explicitQI = explicitQ_.boundaryField()[patchI];
    const vectorField& explicitMI = explicitM_.boundaryField()[patchI];
    const vectorField& explicitMQI = explicitMQ_.boundaryField()[patchI];

    const tensorField& RMI = RM_.internalField();
    const tensorField& refLambdaI = refLambda_.internalField();
    const unallocLabelList& fc = mesh().boundary()[patchI].faceCells();
    
    forAll(cu, faceI)
    {
        ////// W equation

        // W part (Laplacian)
        // upper
        cu[faceI](0,0) += CQWI[faceI].xx()/deltaf[faceI];
        cu[faceI](0,1) += CQWI[faceI].xy()/deltaf[faceI];
        cu[faceI](0,2) += CQWI[faceI].xz()/deltaf[faceI];

        cu[faceI](1,0) += CQWI[faceI].yx()/deltaf[faceI];
        cu[faceI](1,1) += CQWI[faceI].yy()/deltaf[faceI];
        cu[faceI](1,2) += CQWI[faceI].yz()/deltaf[faceI];

        cu[faceI](2,0) += CQWI[faceI].zx()/deltaf[faceI];
        cu[faceI](2,1) += CQWI[faceI].zy()/deltaf[faceI];
        cu[faceI](2,2) += CQWI[faceI].zz()/deltaf[faceI];

        // lower
        cl[faceI](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        cl[faceI](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        cl[faceI](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        cl[faceI](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        cl[faceI](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        cl[faceI](1,2) += -CQWI[faceI].yz()/deltaf[faceI];

        cl[faceI](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        cl[faceI](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        cl[faceI](2,2) += -CQWI[faceI].zz()/deltaf[faceI];

        // diag
        d[fc[faceI]](0,0) += -CQWI[faceI].xx()/deltaf[faceI];
        d[fc[faceI]](0,1) += -CQWI[faceI].xy()/deltaf[faceI];
        d[fc[faceI]](0,2) += -CQWI[faceI].xz()/deltaf[faceI];

        d[fc[faceI]](1,0) += -CQWI[faceI].yx()/deltaf[faceI];
        d[fc[faceI]](1,1) += -CQWI[faceI].yy()/deltaf[faceI];
        d[fc[faceI]](1,2) += -CQWI[faceI].yz()/deltaf[faceI];

        d[fc[faceI]](2,0) += -CQWI[faceI].zx()/deltaf[faceI];
        d[fc[faceI]](2,1) += -CQWI[faceI].zy()/deltaf[faceI];
        d[fc[faceI]](2,2) += -CQWI[faceI].zz()/deltaf[faceI];

        
        // Theta part
        
        // upper
        cu[faceI](0,3) += (1-wf[faceI])*CQThetaI[faceI].xx();
        cu[faceI](0,4) += (1-wf[faceI])*CQThetaI[faceI].xy();
        cu[faceI](0,5) += (1-wf[faceI])*CQThetaI[faceI].xz();

        cu[faceI](1,3) += (1-wf[faceI])*CQThetaI[faceI].yx();
        cu[faceI](1,4) += (1-wf[faceI])*CQThetaI[faceI].yy();
        cu[faceI](1,5) += (1-wf[faceI])*CQThetaI[faceI].yz();
        
        cu[faceI](2,3) += (1-wf[faceI])*CQThetaI[faceI].zx();
        cu[faceI](2,4) += (1-wf[faceI])*CQThetaI[faceI].zy();
        cu[faceI](2,5) += (1-wf[faceI])*CQThetaI[faceI].zz();

        // lower
        cl[faceI](0,3) += wf[faceI]*CQThetaI[faceI].xx();
        cl[faceI](0,4) += wf[faceI]*CQThetaI[faceI].xy();
        cl[faceI](0,5) += wf[faceI]*CQThetaI[faceI].xz();
        
        cl[faceI](1,3) += wf[faceI]*CQThetaI[faceI].yx();
        cl[faceI](1,4) += wf[faceI]*CQThetaI[faceI].yy();
        cl[faceI](1,5) += wf[faceI]*CQThetaI[faceI].yz();
        
        cl[faceI](2,3) += wf[faceI]*CQThetaI[faceI].zx();
        cl[faceI](2,4) += wf[faceI]*CQThetaI[faceI].zy();
        cl[faceI](2,5) += wf[faceI]*CQThetaI[faceI].zz();

        // diag
        d[fc[faceI]](0,3) += wf[faceI]*CQThetaI[faceI].xx();
        d[fc[faceI]](0,4) += wf[faceI]*CQThetaI[faceI].xy();
        d[fc[faceI]](0,5) += wf[faceI]*CQThetaI[faceI].xz();
        
        d[fc[faceI]](1,3) += wf[faceI]*CQThetaI[faceI].yx();
        d[fc[faceI]](1,4) += wf[faceI]*CQThetaI[faceI].yy();
        d[fc[faceI]](1,5) += wf[faceI]*CQThetaI[faceI].yz();
        
        d[fc[faceI]](2,3) += wf[faceI]*CQThetaI[faceI].zx();
        d[fc[faceI]](2,4) += wf[faceI]*CQThetaI[faceI].zy();
        d[fc[faceI]](2,5) += wf[faceI]*CQThetaI[faceI].zz();

        // source
        source[fc[faceI]](0) -= explicitQI[faceI].x();
        source[fc[faceI]](1) -= explicitQI[faceI].y();
        source[fc[faceI]](2) -= explicitQI[faceI].z();


        ////// Theta equation

        // Laplacian part

        // upper
        cu[faceI](3,3) += CMThetaI[faceI].xx()/deltaf[faceI];
        cu[faceI](3,4) += CMThetaI[faceI].xy()/deltaf[faceI];
        cu[faceI](3,5) += CMThetaI[faceI].xz()/deltaf[faceI];
        
        cu[faceI](4,3) += CMThetaI[faceI].yx()/deltaf[faceI];
        cu[faceI](4,4) += CMThetaI[faceI].yy()/deltaf[faceI];
        cu[faceI](4,5) += CMThetaI[faceI].yz()/deltaf[faceI];

        cu[faceI](5,3) += CMThetaI[faceI].zx()/deltaf[faceI];
        cu[faceI](5,4) += CMThetaI[faceI].zy()/deltaf[faceI];
        cu[faceI](5,5) += CMThetaI[faceI].zz()/deltaf[faceI];

        // lower
        cl[faceI](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        cl[faceI](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        cl[faceI](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        cl[faceI](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        cl[faceI](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        cl[faceI](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];
        
        cl[faceI](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        cl[faceI](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        cl[faceI](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];

        // diag
        d[fc[faceI]](3,3) += -CMThetaI[faceI].xx()/deltaf[faceI];
        d[fc[faceI]](3,4) += -CMThetaI[faceI].xy()/deltaf[faceI];
        d[fc[faceI]](3,5) += -CMThetaI[faceI].xz()/deltaf[faceI];

        d[fc[faceI]](4,3) += -CMThetaI[faceI].yx()/deltaf[faceI];
        d[fc[faceI]](4,4) += -CMThetaI[faceI].yy()/deltaf[faceI];
        d[fc[faceI]](4,5) += -CMThetaI[faceI].yz()/deltaf[faceI];
        
        d[fc[faceI]](5,3) += -CMThetaI[faceI].zx()/deltaf[faceI];
        d[fc[faceI]](5,4) += -CMThetaI[faceI].zy()/deltaf[faceI];
        d[fc[faceI]](5,5) += -CMThetaI[faceI].zz()/deltaf[faceI];

        // Theta part

        // upper
        cu[faceI](3,3) += (1-wf[faceI])*CMTheta2I[faceI].xx();
        cu[faceI](3,4) += (1-wf[faceI])*CMTheta2I[faceI].xy();
        cu[faceI](3,5) += (1-wf[faceI])*CMTheta2I[faceI].xz();

        cu[faceI](4,3) += (1-wf[faceI])*CMTheta2I[faceI].yx();
        cu[faceI](4,4) += (1-wf[faceI])*CMTheta2I[faceI].yy();
        cu[faceI](4,5) += (1-wf[faceI])*CMTheta2I[faceI].yz();

        cu[faceI](5,3) += (1-wf[faceI])*CMTheta2I[faceI].zx();
        cu[faceI](5,4) += (1-wf[faceI])*CMTheta2I[faceI].zy();
        cu[faceI](5,5) += (1-wf[faceI])*CMTheta2I[faceI].zz();

        // lower
        cl[faceI](3,3) += wf[faceI]*CMTheta2I[faceI].xx();
        cl[faceI](3,4) += wf[faceI]*CMTheta2I[faceI].xy();
        cl[faceI](3,5) += wf[faceI]*CMTheta2I[faceI].xz();

        cl[faceI](4,3) += wf[faceI]*CMTheta2I[faceI].yx();
        cl[faceI](4,4) += wf[faceI]*CMTheta2I[faceI].yy();
        cl[faceI](4,5) += wf[faceI]*CMTheta2I[faceI].yz();

        cl[faceI](5,3) += wf[faceI]*CMTheta2I[faceI].zx();
        cl[faceI](5,4) += wf[faceI]*CMTheta2I[faceI].zy();
        cl[faceI](5,5) += wf[faceI]*CMTheta2I[faceI].zz();

        // diag
        d[fc[faceI]](3,3) += wf[faceI]*CMTheta2I[faceI].xx();
        d[fc[faceI]](3,4) += wf[faceI]*CMTheta2I[faceI].xy();
        d[fc[faceI]](3,5) += wf[faceI]*CMTheta2I[faceI].xz();

        d[fc[faceI]](4,3) += wf[faceI]*CMTheta2I[faceI].yx();
        d[fc[faceI]](4,4) += wf[faceI]*CMTheta2I[faceI].yy();
        d[fc[faceI]](4,5) += wf[faceI]*CMTheta2I[faceI].yz();

        d[fc[faceI]](5,3) += wf[faceI]*CMTheta2I[faceI].zx();
        d[fc[faceI]](5,4) += wf[faceI]*CMTheta2I[faceI].zy();
        d[fc[faceI]](5,5) += wf[faceI]*CMTheta2I[faceI].zz();


        // Explicit part        
        source[fc[faceI]](3) -= explicitMI[faceI].x();
        source[fc[faceI]](4) -= explicitMI[faceI].y();
        source[fc[faceI]](5) -= explicitMI[faceI].z();


        //---- (dr x  Q) term

        // W part

        // upper
        cu[faceI](3,0) += CMQWI[faceI].xx()/deltaf[faceI];
        cu[faceI](3,1) += CMQWI[faceI].xy()/deltaf[faceI];
        cu[faceI](3,2) += CMQWI[faceI].xz()/deltaf[faceI];

        cu[faceI](4,0) += CMQWI[faceI].yx()/deltaf[faceI];
        cu[faceI](4,1) += CMQWI[faceI].yy()/deltaf[faceI];
        cu[faceI](4,2) += CMQWI[faceI].yz()/deltaf[faceI];

        cu[faceI](5,0) += CMQWI[faceI].zx()/deltaf[faceI];
        cu[faceI](5,1) += CMQWI[faceI].zy()/deltaf[faceI];
        cu[faceI](5,2) += CMQWI[faceI].zz()/deltaf[faceI];

        // lower
        cl[faceI](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        cl[faceI](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        cl[faceI](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        cl[faceI](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        cl[faceI](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        cl[faceI](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        cl[faceI](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        cl[faceI](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        cl[faceI](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];

        // diag
        d[fc[faceI]](3,0) += -CMQWI[faceI].xx()/deltaf[faceI];
        d[fc[faceI]](3,1) += -CMQWI[faceI].xy()/deltaf[faceI];
        d[fc[faceI]](3,2) += -CMQWI[faceI].xz()/deltaf[faceI];

        d[fc[faceI]](4,0) += -CMQWI[faceI].yx()/deltaf[faceI];
        d[fc[faceI]](4,1) += -CMQWI[faceI].yy()/deltaf[faceI];
        d[fc[faceI]](4,2) += -CMQWI[faceI].yz()/deltaf[faceI];

        d[fc[faceI]](5,0) += -CMQWI[faceI].zx()/deltaf[faceI];
        d[fc[faceI]](5,1) += -CMQWI[faceI].zy()/deltaf[faceI];
        d[fc[faceI]](5,2) += -CMQWI[faceI].zz()/deltaf[faceI];

        
        // Theta part

        // upper
        cu[faceI](3,3) += (1-wf[faceI])*CMQThetaI[faceI].xx();
        cu[faceI](3,4) += (1-wf[faceI])*CMQThetaI[faceI].xy();
        cu[faceI](3,5) += (1-wf[faceI])*CMQThetaI[faceI].xz();

        cu[faceI](4,3) += (1-wf[faceI])*CMQThetaI[faceI].yx();
        cu[faceI](4,4) += (1-wf[faceI])*CMQThetaI[faceI].yy();
        cu[faceI](4,5) += (1-wf[faceI])*CMQThetaI[faceI].yz();
        
        cu[faceI](5,3) += (1-wf[faceI])*CMQThetaI[faceI].zx();
        cu[faceI](5,4) += (1-wf[faceI])*CMQThetaI[faceI].zy();
        cu[faceI](5,5) += (1-wf[faceI])*CMQThetaI[faceI].zz();

        // lower
        cl[faceI](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        cl[faceI](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        cl[faceI](3,5) += wf[faceI]*CMQThetaI[faceI].xz();
        
        cl[faceI](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        cl[faceI](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        cl[faceI](4,5) += wf[faceI]*CMQThetaI[faceI].yz();
        
        cl[faceI](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        cl[faceI](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        cl[faceI](5,5) += wf[faceI]*CMQThetaI[faceI].zz();

        // diag
        d[fc[faceI]](3,3) += wf[faceI]*CMQThetaI[faceI].xx();
        d[fc[faceI]](3,4) += wf[faceI]*CMQThetaI[faceI].xy();
        d[fc[faceI]](3,5) += wf[faceI]*CMQThetaI[faceI].xz();

        d[fc[faceI]](4,3) += wf[faceI]*CMQThetaI[faceI].yx();
        d[fc[faceI]](4,4) += wf[faceI]*CMQThetaI[faceI].yy();
        d[fc[faceI]](4,5) += wf[faceI]*CMQThetaI[faceI].yz();

        d[fc[faceI]](5,3) += wf[faceI]*CMQThetaI[faceI].zx();
        d[fc[faceI]](5,4) += wf[faceI]*CMQThetaI[faceI].zy();
        d[fc[faceI]](5,5) += wf[faceI]*CMQThetaI[faceI].zz();

        // Explicit part
        vector correctedOwnExplicitMQ = explicitMQI[faceI];

        correctedOwnExplicitMQ =
        (
            refLambdaI[fc[faceI]].T()
          & (RMI[fc[faceI]].T() & correctedOwnExplicitMQ)
        );
        correctedOwnExplicitMQ.x() = 0;
        correctedOwnExplicitMQ =
            (RMI[fc[faceI]] & (refLambdaI[fc[faceI]] & correctedOwnExplicitMQ));

        source[fc[faceI]](3) -= correctedOwnExplicitMQ.x();
        source[fc[faceI]](4) -= correctedOwnExplicitMQ.y();
        source[fc[faceI]](5) -= correctedOwnExplicitMQ.z();

        // Change sign
        cu[faceI] *= -1;
        cl[faceI] *= -1;
    }

    // Info <<  WThetaEqn.coupleLower()[patchI].asSquare() << endl;
    // Info <<  WThetaEqn.coupleUpper()[patchI].asSquare() << endl;
    
    // Pout << "cu: " << cu << endl;
}
