/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "coupledUpdatedLagBeam.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "cubicSpline.H"
#include "spinTensor.H"
#include "permutationTensor.H"
#include "momentBeamRotationFvPatchVectorField.H"
#include "forceBeamDisplacementFvPatchVectorField.H"
#include "axialForceTransverseDisplacementFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledUpdatedLagBeam, 0);
addToRunTimeSelectionTable(beamModel, coupledUpdatedLagBeam, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledUpdatedLagBeam::coupledUpdatedLagBeam
(
    Time& runTime,
    const word& region
)
:
    beamModel(typeName, runTime, region),
    DW_
    (
        IOobject
        (
            "DW",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    W_
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    DTheta_
    (
        IOobject
        (
            "DTheta",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    theta_
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    m_
    (
        IOobject
        (
            "m",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    Q_
    (
        IOobject
        (
            "Q",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("Q", dimForce, vector::zero)
    ),
    dQ_
    (
        IOobject
        (
            "dQ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("dQ", dimForce, vector::zero)
    ),
    M_
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("M", dimForce*dimLength, vector::zero)
    ),
    dM_
    (
        IOobject
        (
            "dM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("dM", dimForce*dimLength, vector::zero)
    ),
    DLambda_
    (
        IOobject
        (
            "DLambda",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    Lambda_
    (
        IOobject
        (
            "Lambda",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    stretchRatio_
    (
        IOobject
        (
            "stretchRatio",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("1", dimless, 1.0)
    ),
    CQDW_
    (
        IOobject
        (
            "CQDW",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea, tensor::zero)
    ),
    CQDTheta_
    (
        IOobject
        (
            "CQDTheta",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea, tensor::zero)
    ),
    CMDTheta_
    (
        IOobject
        (
            "CMDTheta",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea*dimArea, tensor::zero)
    ),
    CMDTheta2_
    (
        IOobject
        (
            "CMDTheta2",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", M_.dimensions(), tensor::zero)
    ),
    dRuDs_
    (
        IOobject
        (
            "dRuDs",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    pMesh_(mesh()),
    pointDW_
    (
        IOobject
        (
            "pointDW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    DWDTheta_
    (
        IOobject
        (
            "DWDTheta",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector6("zero", dimless, vector6::zero)
    ),
    E_(beamProperties().lookup("E")),
    G_(beamProperties().lookup("G")),
    A_("A", dimArea, M_PI*sqr(R().value())),
    I_("I", dimArea*dimArea, M_PI*pow(R().value(), 4)/4),
    J_("J", dimArea*dimArea, M_PI*pow(R().value(), 4)/2),
    EI_(E_*I_),
    GJ_(G_*J_),
    EA_(E_*A_),
    GA_(G_*A_),
    CQ_
    (
        "CQ",
        EA_.dimensions(),
        tensor
        (
            EA_.value(), 0,           0,
            0,           GA_.value(), 0,
            0,           0,           GA_.value()
        )
    ),
    CM_
    (
        "CM",
        EI_.dimensions(),
        tensor
        (
            GJ_.value(), 0,           0,
            0,           EI_.value(), 0,
            0,           0,           EI_.value()
        )
    )
{
    Info << "I: " << I_ << endl;
    Info << "J: " << J_ << endl;
    Info << "EI: " << EI_ << endl;
    Info << "EA: " << EA_ << endl;
    Info << "GA: " << GA_ << endl;

    epsilon_.oldTime();
    Lambda_.oldTime();

    Q_.oldTime();
    M_.oldTime();

#   include "correctCurvedBeam.H"

    volVectorField Ru
    (
        IOobject
        (
            "Ru",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("Ru", dimLength, vector::zero)
    );
    Ru = mesh().C();
    dRuDs_ = fvc::snGrad(Ru);
    dRuDs_ /= mag(dRuDs_);

    // dRuDs_ = vector(1, 0, 0);
    // dRuDs_.boundaryField()[0] *= -1;
    
    // dRuDs_.write();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coupledUpdatedLagBeam::evolve()
{
    Info<< "\nEvolving beam solver for " << mesh().name() << endl;

    const int nCorr
    (
        beamProperties().lookupOrDefault<int>("nCorrectors", 1000)
    );

    const scalar convergenceTol
    (
        beamProperties().lookupOrDefault<scalar>("convergenceTol", 1e-6)
    );

    label iCorr = 0;
    scalar initialResidual = 0;
    scalar currentResidual = 0;
    blockLduMatrix::debug = 0;

    do
    {
        #include "coupledDWDThetaEqn.H"
    }
    while
    (
        (++iCorr < nCorr)
     && (currentResidual > convergenceTol)
    );

    Info << "Initial residual: " << initialResidual
         << ", current residual: " << currentResidual
         << ", iCorr = " << iCorr << endl;

    return initialResidual;
}

void coupledUpdatedLagBeam::updateTotalFields()
{
    W_ += DW_;
    Lambda_ = (DLambda_ & Lambda_.oldTime());

    const surfaceVectorField DWf = fvc::interpolate(DW_);
    const vectorField& DWfI = DWf.internalField();    
    
    const tensorField& DLambdaI = DLambda_.internalField();
    const faceList& faces = mesh().faces();

    const vectorField& points = mesh().points();

    vectorField& pointDWI = pointDW_.internalField();
    
    forAll(DWfI, faceI)
    {
        const face& curFace = faces[faceI];

        vector C0 = curFace.centre(points);

        forAll(curFace, pointI)
        {
            label curPoint = curFace[pointI];

            vector oldR = points[curPoint] - C0;

            vector newR = C0 + DWfI[faceI]
              + (DLambdaI[faceI] & oldR);

            pointDWI[curPoint] = newR - points[curPoint];
        }
    }

    forAll(DWf.boundaryField(), patchI)
    {
        const vectorField& pDWf =
            DWf.boundaryField()[patchI];

        const tensorField& pDLambda =
            DLambda_.boundaryField()[patchI];

        // const tensorField& pInitialLambda =
        //     initialLambda_.boundaryField()[patchI];
        
        const label start =
            mesh().boundaryMesh()[patchI].start();

        forAll(pDWf, faceI)
        {
            const face& curFace = faces[start + faceI];
            vector C0 = curFace.centre(points);

            forAll(curFace, pointI)
            {
                label curPoint = curFace[pointI];

                vector oldR = points[curPoint] - C0;

                vector newR = C0 + pDWf[faceI]
                  + (pDLambda[faceI] & oldR);

                pointDWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    vectorField newPoints = points + pointDWI;

    const_cast<dynamicFvMesh&>(this->mesh()).movePoints(newPoints);

#   include "correctCurvedBeam.H"
    
    volVectorField Ru
    (
        IOobject
        (
            "Ru",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("Ru", dimLength, vector::zero)
    );
    Ru = mesh().C();

    dRuDs_ = fvc::snGrad(Ru);
    dRuDs_ /= mag(dRuDs_);

    // dRuDs_ = vector(1, 0, 0);
}

tmp<vectorField> coupledUpdatedLagBeam::currentBeamPoints(const label bI) const
{
    label nPoints = this->mesh().nCells()+1;

    tmp<vectorField> tCurrentPoints
    (
        new vectorField(nPoints, vector::zero)
    );
    vectorField& curPoints = tCurrentPoints();

    const surfaceVectorField DWf = fvc::interpolate(DW_);
    surfaceVectorField curCf = mesh().Cf() + DWf;
    const vectorField& curCfI = curCf.internalField();

    curPoints[0] = curCf.boundaryField()[startPatchIndex()][0];
    curPoints[nPoints-1] = curCf.boundaryField()[endPatchIndex()][0];
    for (label i=0; i<curCfI.size(); i++)
    {
        curPoints[i+1] = curCfI[i];
    }

    return tCurrentPoints;
}

tmp<vectorField> coupledUpdatedLagBeam::currentDisplacementIncrement() const
{
    label nPoints = this->mesh().nCells() + 1;

    tmp<vectorField> tDW
    (
        new vectorField(nPoints, vector::zero)
    );
    vectorField& DW = tDW();

    const surfaceVectorField DWf = fvc::interpolate(DW_);

    const vectorField& DWfI = DWf.internalField();

    DW[0] = DWf.boundaryField()[startPatchIndex()][0];
    DW[nPoints-1] = DWf.boundaryField()[endPatchIndex()][0];
    for (label i=0; i<DWfI.size(); i++)
    {
        DW[i+1] = DWfI[i];
    }

    return tDW;
}

tmp<tensorField> coupledUpdatedLagBeam::currentRotationIncrement() const
{
    label nPoints = this->mesh().nCells() + 1;

    tmp<tensorField> tDLambda
    (
        new tensorField(nPoints, tensor::zero)
    );
    tensorField& DLambda = tDLambda();

    const surfaceTensorField& DLambdaf = DLambda_;
    const tensorField& DLambdafI = DLambdaf.internalField();

    DLambda[0] = DLambdaf.boundaryField()[startPatchIndex()][0];
    DLambda[nPoints-1] = DLambdaf.boundaryField()[endPatchIndex()][0];
    for (label i=0; i<DLambdafI.size(); i++)
    {
        DLambda[i+1] = DLambdafI[i];
    }

    return tDLambda;
}
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
