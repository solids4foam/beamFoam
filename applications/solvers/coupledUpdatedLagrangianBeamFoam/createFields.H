    Info<< "Reading field DW\n" << endl;

    volVectorField DW
    (
        IOobject
        (
            "DW",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volVectorField W
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("W", dimLength, vector::zero)
    );

    Info<< "Reading field DTheta\n" << endl;

    volVectorField DTheta
    (
        IOobject
        (
            "DTheta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field q\n" << endl;

    volVectorField q
    (
        IOobject
        (
            "q",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field m\n" << endl;

    volVectorField m
    (
        IOobject
        (
            "m",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading beam properties\n" << endl;

    IOdictionary beamProperties
    (
        IOobject
        (
            "beamProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading beam cross-section radius\n" << endl;

    dimensionedScalar R
    (
        beamProperties.lookup("R")
    );

    Info<< "Reading beam Young modulus\n" << endl;

    dimensionedScalar E
    (
        beamProperties.lookup("E")
    );

    Info<< "Reading beam shear modulus\n" << endl;

    dimensionedScalar G
    (
        beamProperties.lookup("G")
    );

    dimensionedScalar A("A", dimArea, M_PI*sqr(R.value()));

    dimensionedScalar I("I", dimArea*dimArea, M_PI*pow(R.value(), 4)/4);
    dimensionedScalar J("J", dimArea*dimArea, M_PI*pow(R.value(), 4)/2);

    dimensionedScalar EI = E*I;
    dimensionedScalar GJ = G*J;
    dimensionedScalar EA = E*A;
    dimensionedScalar GA = G*A;

    Info << "I: " << I << endl;
    Info << "J: " << J << endl;
    Info << "EI: " << EI << endl;
    Info << "EA: " << EA << endl;
    Info << "GA: " << GA << endl;


    label startPatchIndex = -1;
    label endPatchIndex = -1;
    {
        word startPatchName(beamProperties.lookup("startPatchName"));

        polyPatchID startPatch(startPatchName, mesh.boundaryMesh());

        if (!startPatch.active())
        {
            FatalErrorIn(args.executable())
              << "Patch name " << startPatchName << " not found."
              << abort(FatalError);
        }

        startPatchIndex = startPatch.index();
        
        word endPatchName(beamProperties.lookup("endPatchName"));

        polyPatchID endPatch(endPatchName, mesh.boundaryMesh());

        if (!endPatch.active())
        {
            FatalErrorIn(args.executable())
              << "Patch name " << endPatchName << " not found."
              << abort(FatalError);
        }

        endPatchIndex = endPatch.index();
    }

    // tensor CQ
    // (
    //     EA.value(), 0,          0,
    //     0,          GA.value(), 0,
    //     0,          0,          GA.value()
    // );

    // tensor CTheta
    // (
    //     0, 0, 0,
    //     0, 0,-1,
    //     0, 1, 0
    // );

dimensionedTensor CQ
(
    "CQ",
    EA.dimensions(),
    tensor
    (
        EA.value(), 0,          0,
        0,          GA.value(), 0,
        0,          0,          GA.value()
    )
);

Info << "CQ: " << CQ << endl;

// tensor CM
// (
//     GJ.value(), 0,          0,
//     0,          EI.value(), 0,
//     0,          0,          EI.value()
// );

dimensionedTensor CM
(
    "CM",
    EI.dimensions(),
    tensor
    (
        GJ.value(), 0,          0,
        0,          EI.value(), 0,
        0,          0,          EI.value()
    )
);

Info << "CM: " << CM << endl;

dimensionedVector t ("t", dimless, vector(1,0,0));
surfaceScalarField magSf = (t & mesh.Sf());
surfaceVectorField N = mesh.Sf()/mesh.magSf();

surfaceScalarField sign = magSf;
sign /= mag(sign);

// Allow solution in all three direction
Vector<label>& solutionD = const_cast<Vector<label>&>(mesh.solutionD());
solutionD = Vector<label>(1, 1, 1);

dimensionedTensor II("II", dimless, tensor::I);
dimensionedVector ii("ii", dimless, vector(1, 0, 0));

surfaceVectorField Q
(
    IOobject
    (
        "Q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    CQ & fvc::snGrad(DW)
);
Q.oldTime();

surfaceVectorField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("dQ", Q.dimensions(), vector::zero)
);

surfaceVectorField M
(
    IOobject
    (
        "M",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    CM & fvc::snGrad(DTheta)
);
M.oldTime();

surfaceVectorField dM
(
    IOobject
    (
        "dM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("dM", M.dimensions(), vector::zero)
);

// dimensionedTensor IIw("IIw", Cw.dimensions(), tensor::I);
// dimensionedTensor IItheta("IItheta", Ctheta.dimensions(), tensor::I);

surfaceTensorField DLambda
(
    IOobject
    (
        "DLambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedTensor("I", dimless, tensor::I)
);

surfaceTensorField CQDW
(
    IOobject
    (
        "CQDW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    CQ
    // dimensionedTensor("CQDW", EA.dimensions(), CQ)
);

surfaceTensorField CMDTheta
(
    IOobject
    (
        "CMDTheta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    CM
    // dimensionedTensor("CMDTheta", EI.dimensions(), CM)
);

surfaceTensorField CMDTheta2
(
    IOobject
    (
        "CMDTheta2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   -mulByPermutationTensor(M.oldTime())
);

volVectorField Ru
(
    IOobject
    (
        "Ru",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("R0", dimLength, vector::zero)
);

Ru = mesh.C();

surfaceVectorField dRuDs
(
    IOobject
    (
        "dRuDs",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::snGrad(Ru)
);
dRuDs.write();

#include "correctCurvedBeam.H"

surfaceTensorField CQDTheta
(
    IOobject
    (
        "CQDTheta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
   -(CQDW & mulByPermutationTensor(dRuDs))
);

pointMesh pMesh(mesh);

pointVectorField pointDW
(
    IOobject
    (
        "pointDW",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("0", dimLength, vector::zero)
);

// Block vector field for displacement (first entry) and
// rotation (second entry).
Info << "Creating field DWDTheta\n" << endl;
volVector6Field DWDTheta
(
    IOobject
    (
        "DWDTheta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector6("zero", dimless, vector6::zero)
);
