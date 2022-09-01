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

#include "coupledTotalLagBeam.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "cubicSpline.H"
#include "spinTensor.H"
#include "pseudoVector.H"
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

defineTypeNameAndDebug(coupledTotalLagBeam, 0);
addToRunTimeSelectionTable(beamModel, coupledTotalLagBeam, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledTotalLagBeam::coupledTotalLagBeam
(
    Time& runTime,
    const word& region
)
:
    beamModel(typeName, runTime, region),
    W_
    (
        IOobject
        (
            "W",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    Theta_
    (
        IOobject
        (
            "Theta",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh()
    ),
    refTheta_
    (
        IOobject
        (
            "refTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    // m_
    // (
    //     IOobject
    //     (
    //         "m",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh()
    // ),
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
        dimensionedVector("0", dimForce, vector::zero)
    ),
    explicitQ_
    (
        IOobject
        (
            "explicitQ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce, vector::zero)
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
    explicitM_
    (
        IOobject
        (
            "explicitM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimForce*dimLength, vector::zero)
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
    refLambda_
    (
        IOobject
        (
            "refLambda",
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
    CQW_
    (
        IOobject
        (
            "CQW",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea, tensor::zero)
    ),
    CQTheta_
    (
        IOobject
        (
            "CQTheta",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea, tensor::zero)
    ),
    CMTheta_
    (
        IOobject
        (
            "CMTheta",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea*dimArea, tensor::zero)
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
    i_
    (
        IOobject
        (
            "i",
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
    pointW_
    (
        IOobject
        (
            "pointW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh_,
        dimensionedVector("0", dimLength, vector::zero)
    ),
    WTheta_
    (
        IOobject
        (
            "WTheta",
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
    
    // Calc element lengths
    {
        const fvMesh& mesh = this->mesh();
        label nCellZones = mesh.cellZones().size();

        if (nCellZones < 2)
        {
            const surfaceVectorField& Cf = mesh.Cf();
            const vectorField& CfI = Cf.internalField();

            label nPoints = mesh.nCells()+1;

            vectorField beamPoints(nPoints, vector::zero);
            beamPoints[0] = Cf.boundaryField()[startPatchIndex()][0];
            beamPoints[nPoints-1] = Cf.boundaryField()[endPatchIndex()][0];
            for (label i=0; i<CfI.size(); i++)
            {
                beamPoints[i+1] = Cf[i];
            }

            cubicSpline spline
            (
                beamPoints,
                cubicSpline::CLAMPED_CALC,
                vector(0, 0, 0),
                cubicSpline::CLAMPED_CALC,
                vector(0, 0, 0)
            );

            // Set segment lengths
            this->L().internalField() = spline.segLengths();
        }
        else
        {
            for (label i=0; i<nCellZones; i++)
            {
                vectorField curBeamPoints = points(i);

                cubicSpline spline
                (
                    curBeamPoints,
                    cubicSpline::CLAMPED_CALC,
                    vector(0, 0, 0),
                    cubicSpline::CLAMPED_CALC,
                    vector(0, 0, 0)
                );

                const labelList& curBeamCells = mesh.cellZones()[i];
                scalarField curSegLengths = spline.segLengths();

                forAll(curBeamCells, cellI)
                {
                    // Set segment lengths
                    label curCell = curBeamCells[cellI];
                    this->L().internalField()[curCell] = curSegLengths[cellI];
                }
            }
        }
    }

    // It is assumed that the beam is aligned with x-axis
    // This is valid even for initialy curved beams duo to
    // transformation which will be applied
    i_ = vector(1, 0, 0);
    i_.boundaryField()[startPatchIndex()] *= -1;
    dRuDs_ = vector(1, 0, 0);
    dRuDs_.boundaryField()[startPatchIndex()] *= -1;

    // dRuDs_.write();

    if (contactActive())
    {
        W_.storePrevIter();
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coupledTotalLagBeam::evolve()
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
    scalar curConvergenceTol = convergenceTol;

    const scalar relConvergenceTol
    (
        beamProperties().lookupOrDefault<scalar>("relConvergenceTol", 0)
    );

    const bool debug
    (
        beamProperties().lookupOrDefault<bool>("debug", false)
    );

    scalar initialResidual = 0;
    scalar currentResidual = 0;
    blockLduMatrix::debug = debug;

    label iContactCorr = 0;
    label nContactCorr = 1;
    scalar contactConvTol = 1;
    scalar curContactResidual = 1;
    if (contactActive())
    {
        // nContactCorr = contact().nCorr();
        // contactConvTol = contact().convergenceTol();
    }

    do
    {
        iOuterCorr() = 0;
        do
        {
            #include "coupledWThetaEqn.H"

            // Info << "avg W: "
            //     << max(mag(W_.internalField().component(2))) << endl;

            curConvergenceTol = initialResidual*relConvergenceTol;
            if (curConvergenceTol < convergenceTol)
            {
                curConvergenceTol = convergenceTol;
            }
        }
        while
        (
            (++iOuterCorr() < nCorr)
         && (currentResidual > curConvergenceTol)
        );
        
        Info << "Initial residual: " << initialResidual
             << ", current residual: " << currentResidual
             << ", iCorr = " << iOuterCorr() << endl;

        if (contactActive())
        {
            curContactResidual = contact().update();
            W_.storePrevIter();
            Info << "curContactResidual: "
                << curContactResidual << endl; 
        }
    }
    while
    (
        (++iContactCorr < nContactCorr)
     && (curContactResidual > contactConvTol)
     // && (initialResidual > contactConvTol)
    );

    if (contactActive())
    {
        Info << "Number of contact corrections: "
             << iContactCorr  << endl;
        // contact().updatePenaltyLawOffset();
    }

    return initialResidual;
}

void coupledTotalLagBeam::updateTotalFields()
{
    const surfaceVectorField Wf = fvc::interpolate(W_);
    const vectorField& WfI = Wf.internalField();    

    const tensorField& LambdaI = Lambda_.internalField();
    const tensorField& refLambdaI = refLambda_.internalField();
    const faceList& faces = mesh().faces();

    const vectorField& points = mesh().points();

    vectorField& pointWI = pointW_.internalField();

    forAll(WfI, faceI)
    {
        const face& curFace = faces[faceI];

        vector C0 = curFace.centre(points);

        forAll(curFace, pointI)
        {
            label curPoint = curFace[pointI];

            vector oldR = points[curPoint] - C0;

            vector newR = C0 + WfI[faceI]
              + (LambdaI[faceI] & (refLambdaI[faceI].T() & oldR));

            pointWI[curPoint] = newR - points[curPoint];
        }
    }

    forAll(Wf.boundaryField(), patchI)
    {
        const vectorField& pWf =
            Wf.boundaryField()[patchI];

        const tensorField& pLambda =
            Lambda_.boundaryField()[patchI];

        const tensorField& pRefLambda =
            refLambda_.boundaryField()[patchI];
        
        const label start =
            mesh().boundaryMesh()[patchI].start();

        forAll(pWf, faceI)
        {
            const face& curFace = faces[start + faceI];
            vector C0 = curFace.centre(points);

            forAll(curFace, pointI)
            {
                label curPoint = curFace[pointI];

                vector oldR = points[curPoint] - C0;

                vector newR = C0 + pWf[faceI]
                  + (pLambda[faceI] & (pRefLambda[faceI].T() & oldR));

                pointWI[curPoint] = newR - points[curPoint];
            }
        }
    }
}

tmp<vectorField> coupledTotalLagBeam::currentBeamPoints(const label bI) const
{
    const fvMesh& mesh = this->mesh();
  
    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<vectorField> tCurrentPoints
        (
            new vectorField(nPoints, vector::zero)
        );
        vectorField& curPoints = tCurrentPoints();

        const surfaceVectorField Wf = fvc::interpolate(W_);
        surfaceVectorField curCf = mesh.Cf() + Wf;
        const vectorField& curCfI = curCf.internalField();

        curPoints[0] = curCf.boundaryField()[startPatchIndex()][0];
        curPoints[nPoints-1] = curCf.boundaryField()[endPatchIndex()][0];
        for (label i=0; i<curCfI.size(); i++)
        {
            curPoints[i+1] = curCfI[i];
        }

        return tCurrentPoints;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];

        label nPoints = cz.size() + 1;

        tmp<vectorField> tCurrentPoints
        (
            new vectorField(nPoints, vector::zero)
        );
        vectorField& curPoints = tCurrentPoints();

        const surfaceVectorField Wf = fvc::interpolate(W_);
        surfaceVectorField curCf = mesh.Cf() + Wf;
        const vectorField& curCfI = curCf.internalField();

        const labelList startPatchCells =
            mesh.boundary()[startPatchIndex()].faceCells();

        forAll(curCf.boundaryField()[startPatchIndex()], faceI)
        {
            if (cz.whichCell(startPatchCells[faceI]) != -1)
            {
                curPoints[0] =
                    curCf.boundaryField()[startPatchIndex()][faceI];
            }
        }

        const labelList endPatchCells =
            mesh.boundary()[endPatchIndex()].faceCells();

        forAll(curCf.boundaryField()[startPatchIndex()], faceI)
        {
            if (cz.whichCell(endPatchCells[faceI]) != -1)
            {
                curPoints[nPoints-1] =
                    curCf.boundaryField()[endPatchIndex()][faceI];
            }
        }

        const labelList& nei = mesh.neighbour();

        label index = 1;
        for (label i=0; i<curCfI.size(); i++)
        {
            if (cz.whichCell(nei[i]) != -1)
            {
                curPoints[index] = curCfI[i];
                index++;
            }
        }

        return tCurrentPoints;
    }
}

tmp<vectorField> coupledTotalLagBeam::currentDisplacementIncrement() const
{
    label nPoints = this->mesh().nCells() + 1;

    tmp<vectorField> tDW
    (
        new vectorField(nPoints, vector::zero)
    );
    vectorField& DW = tDW();

    const surfaceVectorField DWf =
        fvc::interpolate(W_)
      - fvc::interpolate(W_.oldTime());

    const vectorField& DWfI = DWf.internalField();

    DW[0] = DWf.boundaryField()[startPatchIndex()][0];
    DW[nPoints-1] = DWf.boundaryField()[endPatchIndex()][0];
    for (label i=0; i<DWfI.size(); i++)
    {
        DW[i+1] = DWfI[i];
    }

    return tDW;
}
  
tmp<tensorField> coupledTotalLagBeam::currentRotationIncrement() const
{
    label nPoints = this->mesh().nCells() + 1;

    tmp<tensorField> tDLambda
    (
        new tensorField(nPoints, tensor::zero)
    );
    tensorField& DLambda = tDLambda();

    const surfaceTensorField DLambdaf = (Lambda_ & inv(Lambda_.oldTime()));

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
