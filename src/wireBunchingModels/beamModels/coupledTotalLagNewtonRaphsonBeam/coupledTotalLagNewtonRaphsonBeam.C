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

#include "coupledTotalLagNewtonRaphsonBeam.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "cubicSpline.H"
#include "HermiteSpline.H"
#include "spinTensor.H"
#include "pseudoVector.H"
#include "permutationTensor.H"
#include "momentBeamRotationFvPatchVectorField.H"
#include "momentBeamRotationNRFvPatchVectorField.H"
#include "forceBeamDisplacementFvPatchVectorField.H"
#include "forceBeamDisplacementNRFvPatchVectorField.H"
#include "followerForceBeamDisplacementNRFvPatchVectorField.H"
#include "axialForceTransverseDisplacementFvPatchVectorField.H"
#include "axialForceTransverseDisplacementNRFvPatchVectorField.H"
#include "extrapolatedBeamRotationFvPatchVectorField.H"

#include "mergePoints.H"
#include "scalarMatrices.H"
#include "denseMatrixHelperFunctions.H"
#include "BlockEigenSolverOF.H"

// #include "beamHelperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(coupledTotalLagNewtonRaphsonBeam, 0);
addToRunTimeSelectionTable
(
    beamModel,
    coupledTotalLagNewtonRaphsonBeam,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledTotalLagNewtonRaphsonBeam::coupledTotalLagNewtonRaphsonBeam
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
    WIncrement_
    (
        IOobject
        (
            "WIncrement",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", W_.dimensions(), vector::zero)
    ),
    U_
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
         //   fvc::ddt(W_)
 	      dimensionedVector("0", W_.dimensions()/dimTime, vector::zero)
      ),
        Accl_
        (
 	      IOobject
 	  (
 	    "Accl",
 	    runTime.timeName(),
 	    mesh(),
 	    IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
 	  ),
 	  mesh(),
 	  dimensionedVector("0", U_.dimensions()/dimTime, vector::zero)
        ),
    DW_
    (
        IOobject
        (
            "DW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", W_.dimensions(), vector::zero)
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
    ThetaIncrement_
    (
        IOobject
        (
            "ThetaIncrement",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", Theta_.dimensions(), vector::zero)
    ),
    Omega_
    (
        IOobject
        (
            "Omega",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
      //   fvc::ddt(Theta_)
 	    dimensionedVector("0", Theta_.dimensions()/dimTime, vector::zero)
      ),
      dotOmega_
      (
         IOobject
         (
             "dotOmega",
             runTime.timeName(),
             mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
 	        mesh(),
         dimensionedVector("0", Omega_.dimensions()/dimTime, vector::zero)
    ),
    DTheta_
    (
        IOobject
        (
            "DTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", Theta_.dimensions(), vector::zero)
    ),
    updatedLagrangian_
    (
        beamProperties().lookupOrDefault<bool>("updatedLagrangian", false)
    ),
    kirchhoffBeam_
    (
        beamProperties().lookupOrDefault<bool>("kirchhoffBeam", false)
    ),
    totW_
    (
        IOobject
        (
            "totW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
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
    Qa_
    (
        IOobject
        (
            "Qa",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("0", dimForce, 0)
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
    Mref_
    (
        IOobject
        (
            "Mref",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("Mref", dimForce*dimLength, vector::zero)
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
    explicitMQ_
    (
        IOobject
        (
            "explicitMQ",
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
    RM_
    (
        IOobject
        (
            "RM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    refWf_
    (
        IOobject
        (
            "refWf",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
    ),
    refW_
    (
        IOobject
        (
            "refW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimLength, vector::zero)
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
    refRM_
    (
        IOobject
        (
            "refRM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, tensor::I)
    ),
    refTangent_
    (
        IOobject
        (
            "refTangent",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("x-axis", dimless, vector(1, 0, 0))
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimForce*dimLength, tensor::zero)
    ),
    CMTheta_
    (
        IOobject
        (
            "CMTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimPressure*dimArea*dimArea, tensor::zero)
    ),
    CMTheta2_
    (
        IOobject
        (
            "CMTheta2",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", M_.dimensions(), tensor::zero)
    ),
    CMQW_
    (
        IOobject
        (
            "CMQW",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimForce*dimLength, tensor::zero)
    ),
    CMQTheta_
    (
        IOobject
        (
            "CMQTheta",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimForce*dimLength, tensor::zero)
    ),
    dR0Ds_
    (
        IOobject
        (
            "dR0Ds",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    Gamma_
    (
        IOobject
        (
            "Gamma",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless, vector::zero)
    ),
    kirchhoffTransTensor_
    (
        IOobject
        (
            "kirchhoffTransTensor",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("0", dimless, tensor::I)
    ),
    K_
    (
        IOobject
        (
            "K",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimless/dimLength, vector::zero)
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
    // WTheta_
    // (
    //     IOobject
    //     (
    //         "WTheta",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedVector6("zero", dimless, vector6(1e-6))
    // ),
    totalIter_(0),
  //  nCV_(beamProperties().lookup("nSegments")),
    // E_(beamProperties().lookup("E")),
    // G_(beamProperties().lookup("G")),
    // rho_("rho", dimDensity, 0),
    // A_("A", dimArea, M_PI*sqr(R().value())),
    // I_("I", dimArea*dimArea, M_PI*pow(R().value(), 4)/4),
    // J_("J", dimArea*dimArea, M_PI*pow(R().value(), 4)/2),
    // EI_(E_*I_),
    // GJ_(G_*J_),
    // EA_(E_*A_),
    // GA_(G_*A_),
    // CQ_
    // (
    //     "CQ",
    //     EA().dimensions(),
    //     tensor
    //     (
    //         EA().value(), 0,           0,
    //         0,           GA().value(), 0,
    //         0,           0,           GA().value()
    //     )
    // ),
    CQ_
    (
        IOobject
        (
            "CQ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        // mesh(),
        indicator(0)
       *dimensionedTensor
        (
            "CQ",
            EA(0).dimensions(),
            tensor
            (
                EA(0).value(), 0, 0,
                0, GA(0).value(), 0,
                0, 0, GA(0).value()
            )
        )
    ),
    CDQ_
    (
        IOobject
        (
            "CDQ",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        // mesh(),
        indicator(0)
       *dimensionedTensor
        (
            "CQ",
            EA().dimensions(),
            tensor
            (
                EA().value(), 0, 0,
                0, GA().value(), 0,
                0, 0, GA().value()
            )
        )
        // CQ_
    ),
    CM_
    (
        IOobject
        (
            "CM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        // mesh(), //CM_
        indicator(0)
       *dimensionedTensor
        (
            "CM",
            EI().dimensions(),
            tensor
            (
                GJ().value(), 0,              0,
                0,            EIyy().value(), 0,
                0,            0,              EIzz().value()
            )
        )
    ),
    // CM_
    // (
    //     "CM",
    //     EI().dimensions(),
    //     diagTensor
    //     (
    //         GJ().value(),
    //         EI().value(),
    //         EI().value()
    //     )
    //     // (
    //     //     GJ().value(), 0,           0,
    //     //     0,           EI().value(), 0,
    //     //     0,           0,           EI().value()
    //     // )
    // ),
    CDM_
    (
        IOobject
        (
            "CDM",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(), //CM_
        dimensionedTensor
        (
            "CDM",
            EI().dimensions(),
            tensor
            (
                GJ().value(), 0,              0,
                0,            EIyy().value(), 0,
                0,            0,              EIzz().value()
            )
        )
    ),
    CDMDGamma_
    (
        IOobject
        (
            "CDMDGamma",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(), //CM_
        dimensionedTensor
        (
            "CDMDGamma",
            dimForce*dimLength,
            tensor::zero
        )
    ),
    CDQDK_
    (
        IOobject
        (
            "CDQDK",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor
        (
            "CDQDK",
            dimForce*dimLength,
            tensor::zero
        )
    ),
    ARho_
     (
 	  "ARho",
 	  A().dimensions()*rho().dimensions(),
 	  scalar(ARho().value())
     ),
     //   CI_
     //   (
         //   "CI",
         //   I().dimensions(),
         //   tensor
         //   (
             //   kCI_*J().value(),  0,           0,
             //   0,            kCI_*Iyy().value(), 0,
             //   0,            0,           kCI_*Izz().value()
         //   )
     //   ),
      CIRho_
      (
         "CIRho",
         I().dimensions()*rho().dimensions(),
         tensor
         (
             kCI()*JRho().value(),  0,           0,
             0,            kCI()*IyyRho().value(), 0,
             0,            0,           kCI()*IzzRho().value()
         )
      ),
     newmark_(beamProperties().lookupOrDefault<bool>("newmark", false)),
     betaN_(beamProperties().lookupOrDefault<scalar>("newmarkBeta", 0.25)),
     gammaN_(beamProperties().lookupOrDefault<scalar>("newmarkGamma", 0.5)),

    // Drag Force related fields
    dragActive_(beamProperties().lookupOrDefault<bool>("dragActive", false)),
    Cdn_(beamProperties().lookupOrDefault<scalar>("Cdn", 1.0)),
    Cdt_(beamProperties().lookupOrDefault<scalar>("Cdt", 1.0)),


    // ground contact fields
    groundContactActive_(beamProperties().lookupOrDefault<bool>("groundContactActive", false)),
    segNo_(readInt(beamProperties().lookup("nSegments"))),
    gDamping_ (readDouble(beamProperties().lookup("gDamping"))),
    gStiffness_ (readDouble(beamProperties().lookup("gStiffness"))),
    groundZ_ (readDouble(beamProperties().lookup("groundZ"))),

    // Plasticity related fields
    // plasticity_(lookupOrDefault<bool>("plasticity", false)),
    // GammaP_
    // (
    //     IOobject
    //     (
    //         "GammaP",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedVector("0", dimless, vector::zero)
    // ),
    // KP_
    // (
    //     IOobject
    //     (
    //         "KP",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedVector("0", dimless/dimLength, vector::zero)
    // ),
    // plasticityStressResultantReturnPtr_(NULL),
    curvature_
    (
        IOobject
        (
            "curvature",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless/dimLength, 1)
    ),
    totalContactTime_(0),
    totalSolutionTime_(),
    proc_
    (
        IOobject
        (
            "proc",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0)
    )
{
    W_.oldTime();
    U_.oldTime();
    Omega_.oldTime();
    RM_.oldTime();

    Gamma_.oldTime();
    K_.oldTime();
    Lambda_.oldTime();

    Q_.oldTime();
    M_.oldTime();

    // Correct properties for multibeam case
    label nBeams = mesh().cellZones().size();
    if (nBeams > 1)
    {
        // First beam is already initialized
        for (label i=1; i<nBeams; i++)
        {
            CQ_ +=
                indicator(i)
               *dimensionedTensor
                (
                    "CQ",
                    EA(i).dimensions(),
                    tensor
                    (
                        EA(i).value(), 0, 0,
                        0, GA(i).value(), 0,
                        0, 0, GA(i).value()
                    )
                );

            CM_ +=
                indicator(i)
               *dimensionedTensor
                (
                    "CM",
                    EI(i).dimensions(),
                    tensor
                    (
                        GJ(i).value(), 0, 0,
                        0, EIyy(i).value(), 0,
                        0, 0, EIzz(i).value()
                    )
                );
        }

        CDQ_ = CQ_;
        CDM_ = CM_;
    }

    // Plasticity
    // if (beamProperties().found("plasticityModel"))
    // {
    //     // Info << "Plasticity model found" << endl;
    //     plasticityStressResultantReturnPtr_ =
    //         plasticityStressResultantReturn::New
    //         (
    //             beamProperties().lookup("plasticityModel"),
    //             *this
    //         );
    // }

    // GammaP_.oldTime();
    // KP_.oldTime();

    // Calculate tangents if it is not read
    IOobject refTangentHeader
    (
        "refTangent",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    if (!refTangentHeader.typeHeaderOk<surfaceVectorField>(true))
    {
        Info<< "Calculating mean line tangents for initial configuration"
            << endl;

        volVectorField R0
        (
            IOobject
            (
                "R0",
                runTime.timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
             ),
            mesh(),
            dimensionedVector("R0", dimLength, vector::zero)
        );
        R0 = mesh().C();
        //R0.boundaryFieldRef().evaluateCoupled();
        dR0Ds_ = fvc::snGrad(R0);
        dR0Ds_ /= mag(dR0Ds_);

        // Info << mesh().C() << endl;
        // Info << mesh().magSf() << endl;

        refTangent_ = dR0Ds_;

        if (nBeams < 2)
        {
            if (refTangent_.boundaryFieldRef()[startPatchIndex()].size())
            {
                refTangent_.boundaryFieldRef()[startPatchIndex()] *= -1;
            }
            else
            {
                // Find left processor patch
                forAll(mesh().boundary(), patchI)
                {
                    if (isA<processorFvPatch>(mesh().boundary()[patchI]))
                    {
                        const labelList& fc =
                            mesh().boundary()[patchI].faceCells();

                        if (fc[0] == 0)
                        {
                            refTangent_.boundaryFieldRef()[patchI] *= -1;
                        }
                    }
                }
            }
        }
        else
        {
            for (label bI=0; bI<nBeams; bI++)
            {
                const cellZone& cz = mesh().cellZones()[bI];
                label firstCell = min(cz);
                // Pout << bI << ", " << firstCell << endl;

                if (refTangent_.boundaryFieldRef()[startPatchIndex(bI)].size())
                {
                    refTangent_.boundaryFieldRef()[startPatchIndex(bI)] *= -1;
                }
                else
                {
                    // Find left processor patch
                    forAll(mesh().boundary(), patchI)
                    {
                        if (isA<processorFvPatch>(mesh().boundary()[patchI]))
                        {
                            const labelList& fc =
                                mesh().boundary()[patchI].faceCells();

                            // if (Pstream::myProcNo() == 1)
                            // {
                            //     Pout << "fc: " << fc << endl;
                            // }

                            forAll(fc, cI)
                            {
                                if (cz.whichCell(fc[cI]) != -1)
                                {
                                    if (fc[cI] == firstCell)
                                    {
                                        // if (Pstream::myProcNo() == 1)
                                        // {
                                        //     Pout << fc[cI] << " == " << firstCell << endl;
                                        //     Pout << refTangent_.boundaryField()[patchI] << endl;
                                        // }

                                        refTangent_.boundaryFieldRef()[patchI][cI]
                                            *= -1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // Info << refTangent_ << endl;
        // Info << dR0Ds_ << endl;

        //
        // if (Pstream::myProcNo() == 1)
        // {
        //     Pout << refTangent_ << endl;
        // }
        // sleep(5);

        // Pout << "Mean line tangent...done" << endl;
        // sleep(5);

    }
    else
    {
        dR0Ds_ = refTangent_;
        for (label pI=0; pI<nBeams; pI++)
        {
            dR0Ds_.boundaryFieldRef()[startPatchIndex(pI)] *= -1;
        }
    }

    // Info << refTangent_ << endl;

    // Calculate cell-centre reference rotation matrix
    // if it is not read
    IOobject refRMheader
    (
        "refRM",
        runTime.timeName(),
        mesh(),
        IOobject::MUST_READ
    );
    if (!refRMheader.typeHeaderOk<volTensorField>(true))
    {
        Info<< "Calculating cell-centre reference rotation matrix" << endl;

        // Calc cell-centre reference rotation matrix
        interpolateRotationMatrix(*this, refLambda_, refRM_);

        if (debug)
        {
            Info<< refRM_ << endl;
        }
    }


    // Calc element lengths
    if (mesh().nCells() > 0)
    {
        // // For correcting cell centres
        // dynamicFvMesh& m = const_cast<dynamicFvMesh&>(this->mesh());
        // m.movePoints(this->mesh().points());
        // m.moving(false);

        // vectorField& cellCentres =
        //    const_cast<vectorField&>(this->mesh().cellCentres());

        if (nBeams < 2)
        {
	    surfaceVectorField curCf(mesh().Cf() + refWf_);
            vectorField beamPoints(this->beamPointData(curCf));
            vectorField beamTangents(this->beamPointData(refTangent_));

            HermiteSpline spline
            (
                beamPoints,
                beamTangents
            );

            // Set segment lengths
            this->L().primitiveFieldRef() = spline.segLengths();

            // cellCentres = spline.midPoints();

            Pout<< "Beam length: " << spline.length() << endl;
        }
        else
        {
	    surfaceVectorField curCf(mesh().Cf() + refWf_);
            for (label i=0; i<nBeams; i++)
            {
                vectorField beamPoints(this->beamPointData(curCf, i));
                vectorField beamTangents(this->beamPointData(refTangent_, i));

                if (beamPoints.size())
                {
                    // Pout << beamPoints << endl;
                    // Pout << beamTangents << endl;
                    // Pout << i << " Length ... 0" << endl;
                    // sleep(5);

                    HermiteSpline spline
                    (
                        beamPoints,
                        beamTangents
                    );

                    // Pout << i << " Length ... 1" << endl;
                    // sleep(5);

                    const labelList& curBeamCells = mesh().cellZones()[i];
                    scalarField curSegLengths(spline.segLengths());
                    // vectorField newCellCentres = spline.midPoints();

                    // Pout << i << " Length ... 2" << endl;
                    // sleep(5);

                    forAll(curBeamCells, cellI)
                    {
                        // Set segment lengths
                        label curCell = curBeamCells[cellI];
                        this->L().primitiveFieldRef()[curCell] =
                            curSegLengths[cellI];

                        // // Correct cell centres
                        // cellCentres[curCell] = newCellCentres[cellI];
                    }
                }
            }
        }
    }
    Pout<< __FILE__ << " " << __LINE__ << endl;

    W_.storePrevIter();
    Theta_.storePrevIter();


    // Update contact for restart
    // if (contactActive())
    // {
    //     if (runTime.timeIndex()>1)
    //     {
	//     // Info << "Inside if" << endl;
    //         contact().update();
    //         contact().finalUpdate();
    //     }
    //     else
    //     {
	//     // Info << "else loop" << endl;
    //         contact();
    //     }
    // }



    // Set Kirchhoff beam transformation tensor
    // if (kirchhoffBeam_)
    // {
    //     // kirchhoffTransTensor_ = vector(1,0,0)*vector(1,0,0);
    //     // Info << kirchhoffTransTensor_ << endl;
    // }


    // Set fields as oriented
    explicitQ_.setOriented(true);
    explicitM_.setOriented(true);
    Q_.setOriented(true);
    M_.setOriented(true);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar coupledTotalLagNewtonRaphsonBeam::evolve()
{
    beamModel::evolve();

    const int nCorr
    (
        beamProperties().lookupOrDefault<int>("nCorrectors", 1000)
    );

    const scalar convergenceTol
    (
        beamProperties().lookupOrDefault<scalar>("convergenceTol", 1e-6)
    );
    scalar curConvergenceTol = convergenceTol;

    const scalar materialTol
    (
        beamProperties().lookupOrDefault<scalar>
        (
            "materialTol",
            curConvergenceTol
        )
    );

    // const scalar relConvergenceTol
    // (
    //     beamProperties().lookupOrDefault<scalar>("relConvergenceTol", 0)
    // );

    const bool debug
    (
        beamProperties().lookupOrDefault<bool>("debug", false)
    );

    scalar initialResidual = 1;
    scalar currentResidual = 1;
    scalar currentMaterialResidual = 0;
    //bool completedElasticPrediction = false;
    //blockLduMatrix::debug = debug;

    scalar curContactResidual = 1;

    iOuterCorr() = 0;
    do
    {
        if (debug)
        {
            Info << "iOuterCorr: " << iOuterCorr() << endl;
        }

        // if (contactActive())
        // {
        //     if (debug)
        //     {
        //         Info << "Updating contact: start" << endl;
        //     }

        //     scalar tStart = runTime().elapsedCpuTime();

        //     // Info << "tstart in CTLNRB file: \n " << tStart << endl;
        //     curContactResidual = contact().update();
        //     scalar tEnd = runTime().elapsedCpuTime();

        //     totalContactTime_ += tEnd - tStart;

        //     if (debug)
        //     {
        //         Pout << "Current total contact update time: "
        //              << totalContactTime_ << endl;
        //     }

        //     if (debug)
        //     {
        //         Info << "curContactResidual: "
        //             << curContactResidual << endl;

        //         Info << "Updating contact: end" << endl;
        //     }
        // }

        scalar tStart = runTime().elapsedCpuTime();
        #include "coupledWThetaEqn_TLNR.H"
        scalar tEnd = runTime().elapsedCpuTime();

        totalSolutionTime_ += tEnd-tStart;

        if (debug)
        {
            Pout<< "Current total solution update time: "
                << totalSolutionTime_ << endl;
        }

        // curConvergenceTol = initialResidual*relConvergenceTol;
        // if (curConvergenceTol < convergenceTol)
        // {
        //     curConvergenceTol = convergenceTol;
        // }
    }
    while
    (
        (++iOuterCorr() < nCorr)
     && (
            (currentResidual > curConvergenceTol)
         || (currentMaterialResidual > materialTol)
        )
    );

    totalIter_ += iOuterCorr();

    Info<< "\nInitial residual: " << initialResidual
        << ", current residual: " << currentResidual
        << ", current material residual: " << currentMaterialResidual
        << ", current contact force residual: " << curContactResidual
        << ",\n iCorr = " << iOuterCorr() << nl
        << "total Iterations " << totalIter_ << endl;

    return initialResidual;
}

void coupledTotalLagNewtonRaphsonBeam::updateTotalFields()
{
    // WIncrement_ = W_ - W_.oldTime(); // updated after each outer iter

    // Calculate rotation angle incrment
    // if (false)
    // {
    //     volTensorField relRM = (RM_ & RM_.oldTime().T());
    //     volVectorField relTheta = pseudoVector(relRM);

    //     // ThetaIncrement_ = relTheta;
    //     ThetaIncrement_ = (RM_.oldTime().T() & relTheta);

    //     Info << "relRM: " << relRM.boundaryField()[1] << endl;

    //     Theta_ = Theta_.oldTime() + ThetaIncrement_;

    //     // ThetaIncrement_ = Theta_ - Theta_.oldTime();

    //     Info << "relTheta[1]: "
    //          << relTheta.boundaryField()[1] << endl;

    //     Info << "ThetaIncrement[1]: "
    //          << ThetaIncrement_.boundaryField()[1] << endl;

    //     Info << "Theta[1]: "
    //          << Theta_.boundaryField()[1] << endl;

    //     Info << "oldTheta[1]: "
    //          << Theta_.oldTime().boundaryField()[1] << endl;

    //     Info << "DTheta[1]: "
    //          << Theta_.boundaryField()[1]
    //           - Theta_.oldTime().boundaryField()[1] << endl;

    //     tensor RM1 = RM_.boundaryField()[1][0];
    //     tensor oldRM1 = RM_.oldTime().boundaryField()[1][0];

    //     tensor RM1star =
    //         rotationMatrix(Theta_.boundaryField()[1][0]);

    //     tensor oldRM1star =
    //         rotationMatrix(Theta_.oldTime().boundaryField()[1][0]);

    //     tensor relRMstar = (RM1star & oldRM1star.T());

    //     Info << "relRMstar: " << relRMstar << endl;

    //     vector relAngle = pseudoVector(relRMstar);
    //     vector AngleIncrement =  (oldRM1star.T() & relAngle);

    //     Info << "RM[1]: " << RM1 << endl;
    //     Info << "RM*[1]: " << RM1star << endl;

    //     Info << "oldRM[1]: " << oldRM1 << endl;
    //     Info << "oldRM*[1]: " << oldRM1star << endl;

    //     Info << "relAngle: " << relAngle << endl;
    //     Info << "AngleIncrement: " << AngleIncrement << endl;
    // }
    // else
    {
        ThetaIncrement_ =  Theta_ - Theta_.oldTime();
    }

    surfaceVectorField Wf(fvc::interpolate(W_) + refWf_);

    //label nCellZones = mesh().cellZones().size();

    const vectorField& WfI = Wf.internalField();

    const tensorField& LambdaI = Lambda_.internalField();
    const tensorField& refLambdaI = refLambda_.internalField();
    const faceList& faces = mesh().faces();

    const vectorField& points = mesh().points();

    vectorField& pointWI(pointW_.primitiveFieldRef());

    forAll(WfI, faceI)
    {
        const face& curFace = faces[faceI];

        vector C0 = curFace.centre(points);

        // Info << faceI << ", " << C0 << endl;

        forAll(curFace, pointI)
        {
            label curPoint = curFace[pointI];

            vector oldR = points[curPoint] - C0;

            // Info << faceI << ", " << oldR << endl;

            vector newR = C0 + WfI[faceI] + //(LambdaI[faceI] & oldR);
                (LambdaI[faceI] & (refLambdaI[faceI] & oldR));

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

                vector newR = C0 + pWf[faceI] +
                 (pLambda[faceI] & (pRefLambda[faceI] & oldR));

                pointWI[curPoint] = newR - points[curPoint];
            }
        }
    }

    // Calc radius of curvature
    // if (false) // Check for parallel runs
    // {
    //     if (nCellZones < 2)
    //     {
    //         vectorField beamPoints = this->currentBeamPoints(0);
    //         vectorField beamTangents = this->currentBeamTangents(0);

    //         HermiteSpline spline
    //         (
    //             beamPoints,
    //             beamTangents
    //         );

    //         const scalarField& segCurvature =
    //             spline.midPointCurvatures();

    //         // Set radius of curvature
    //         curvature_.internalField() = segCurvature;
    //     }
    //     else
    //     {
    //         for (label i=0; i<nCellZones; i++)
    //         {
    //             vectorField beamPoints =
    //                 this->currentBeamPoints(i);
    //             vectorField beamTangents =
    //                 this->currentBeamTangents(i);

    //             HermiteSpline spline
    //             (
    //                 beamPoints,
    //                 beamTangents
    //             );

    //             const labelList& curBeamCells = mesh().cellZones()[i];

    //             const scalarField& segCurvature =
    //                 spline.midPointCurvatures();

    //             forAll(curBeamCells, cellI)
    //             {
    //                 // Set segment lengths
    //                 label curCell = curBeamCells[cellI];
    //                 curvature_.internalField()[curCell] = segCurvature[cellI];
    //             }
    //         }
    //     }
    // }

    // Pout << "updateTotalFields 0" << endl;
    // sleep(5);

    // if (updatedLagrangian_)
    // {
    //     totW_ += W_;
    //     refLambda_ = (Lambda_ & refLambda_);
    //     // refRM_ = (RM_ & refRM_);

    //     // Calc cell-centre reference rotation matrix
    //     interpolateRotationMatrix(*this, refLambda_, refRM_);

    //     volVectorField R0
    //     (
    //         IOobject
    //         (
    //             "R0",
    //             runTime().timeName(),
    //             mesh(),
    //             IOobject::NO_READ,
    //             IOobject::NO_WRITE
    //          ),
    //         mesh(),
    //         dimensionedVector("R0", dimLength, vector::zero)
    //     );
    //     R0 = mesh().C() + W_;
    //     R0.boundaryField().evaluateCoupled();
    //     dR0Ds_ = fvc::snGrad(R0);
    //     dR0Ds_ /= mag(dR0Ds_);
    //     refTangent_ = dR0Ds_;

    //     // Info << "nCellZones: " << nCellZones << endl;
    //     for (label pI=0; pI<nCellZones; pI++)
    //     {
    //         refTangent_.boundaryField()[startPatchIndex(pI)] *= -1;
    //     }

    //     // Info << refTangent_.boundaryField() << endl;

    //     Lambda_ == dimensionedTensor("I", Lambda_.dimensions(), tensor::I);
    //     RM_ == dimensionedTensor("I", RM_.dimensions(), tensor::I);
    //     Gamma_ == dimensionedVector("0", Gamma_.dimensions(), vector::zero);
    //     K_ == dimensionedVector("0", K_.dimensions(), vector::zero);

    //     W_ == dimensionedVector("0", W_.dimensions(), vector::zero);
    //     W_.storePrevIter();

    //     Theta_ == dimensionedVector("0", Theta_.dimensions(), vector::zero);
    //     Theta_.storePrevIter();

    //     // if (contactActive())
    //     // {
    //     //     W_ = dimensionedVector("0", W_.dimensions(), vector::zero);
    //     //     W_.storePrevIter();
    //     // }

    //     vectorField newPoints = points + pointWI;
    //     const_cast<dynamicFvMesh&>(this->mesh()).movePoints(newPoints);

    //     // Calc element lengths and currect cell centres
    //     vectorField& cellCentres =
    //         const_cast<vectorField&>(mesh().cellCentres());
    //     {
    //         label nCellZones = mesh().cellZones().size();

    //         if (nCellZones < 2)
    //         {
    //             vectorField beamPoints = this->beamPointData(mesh().Cf());
    //             vectorField beamTangents = this->beamPointData(refTangent_);

    //             HermiteSpline spline
    //             (
    //                 beamPoints,
    //                 beamTangents
    //             );

    //             // Set segment lengths
    //             this->L().internalField() = spline.segLengths();

    //             cellCentres = spline.midPoints();

    //             Info << "Beam length: " << spline.length() << endl;
    //         }
    //         else
    //         {
    //             for (label i=0; i<nCellZones; i++)
    //             {
    //                 vectorField beamPoints =
    //                     this->beamPointData(mesh().Cf(), i);
    //                 vectorField beamTangents =
    //                     this->beamPointData(refTangent_, i);

    //                 if (beamPoints.size())
    //                 {
    //                     HermiteSpline spline
    //                     (
    //                         beamPoints,
    //                         beamTangents
    //                     );

    //                     const labelList& curBeamCells = mesh().cellZones()[i];
    //                     scalarField curSegLengths = spline.segLengths();
    //                     vectorField newCellCentres = spline.midPoints();

    //                     forAll(curBeamCells, cellI)
    //                     {
    //                         // Set segment lengths
    //                         label curCell = curBeamCells[cellI];
    //                         this->L().internalField()[curCell] =
    //                             curSegLengths[cellI];

    //                         // Correct cell centres
    //                         cellCentres[curCell] = newCellCentres[cellI];
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // Pout << "updateTotalFields 1" << endl;
    // sleep(5);

    // Plasticity
    // if (plasticityStressResultantReturnPtr_.valid())
    // {
    //     plasticityStressResultantReturnPtr_().updateYieldStress();

    //     CDQ_ = CQ_;
    //     CDM_ = CM_;

    //         // dimensionedTensor
    //         // (
    //         //     "CM",
    //         //     CM_.dimensions(),
    //         //     tensor
    //         //     (
    //         //         GJ().value(), 0, 0,
    //         //         0, EI().value(), 0,
    //         //         0, 0, EI().value()
    //         //     )
    //         // );

    //     CDMDGamma_ =
    //         dimensionedTensor("zero", dimForce*dimLength, tensor::zero);

    //     CDQDK_ =
    //         dimensionedTensor("zero", dimForce*dimLength, tensor::zero);

    //     plasticityStressResultantReturnPtr_().writeInfo();
    // }

    // if (conicalPulleys().size())
    // {
    //     label nBeams = contact().splines().size();
    //     label nPulleys = conicalPulleys().size();

    //     for (label bI=0; bI<nBeams; bI++)
    //     {
    //         const scalarPairList& curConicalPulleyContactGaps =
    //             contact().conicalPulleyContactGaps()[bI];

    //         Info << "Beam " << bI << ": "
    //              << curConicalPulleyContactGaps[0];
    //         for (label i=1; i<nPulleys; i++)
    //         {
    //             Info << ", " << curConicalPulleyContactGaps[i];

    //         }
    //         Info << endl;
    //     }
    // }

    // if (toroidalPulleys().size())
    // {
    //     label nBeams = contact().splines().size();
    //     label nPulleys = toroidalPulleys().size();

    //     for (label bI=0; bI<nBeams; bI++)
    //     {
    //         const scalarPairList& curToroidalPulleyContactGaps =
    //             contact().toroidalPulleyContactGaps()[bI];

    //         Info << "Beam " << bI << ": "
    //              << curToroidalPulleyContactGaps[0];
    //         for (label i=1; i<nPulleys; i++)
    //         {
    //             Info << ", " << curToroidalPulleyContactGaps[i];

    //         }
    //         Info << endl;
    //     }
    // }

    beamModel::updateTotalFields();
}


void coupledTotalLagNewtonRaphsonBeam::writeFields()
{
    beamModel::writeFields();

    // Plasticity
    // if (plasticityStressResultantReturnPtr_.valid())
    // {
    //     if (runTime().outputTime())
    //     {
    //         plasticityStressResultantReturnPtr_().writeFields();

    //         // surfaceVectorField dRdS
    //         // (
    //         //     IOobject
    //         //     (
    //         //         "dRdS",
    //         //         runTime.timeName(),
    //         //         mesh(),
    //         //         IOobject::NO_READ,
    //         //         IOobject::AUTO_WRITE
    //         //     ),
    //         //     dR0Ds_ + fvc::snGrad(W_)
    //         // );

    //         // dRdS.write();
    //     }
    // }
}


tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentBeamPoints
(
    const label bI
) const
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

        const surfaceVectorField Wf(fvc::interpolate(W_));
        surfaceVectorField curCf(mesh.Cf() + refWf_ + Wf);

        // Calculate cell-centre mean line tangents
        vectorField tangents(mesh.nCells(), vector::zero);
        const labelList& upperNeiFaces = upperNeiCellFaces();
        const labelList& lowerNeiFaces = lowerNeiCellFaces();
        forAll(tangents, cellI)
        {
            label upperFace = upperNeiFaces[cellI];
            label lowerFace = lowerNeiFaces[cellI];

            if (upperFace > -1)
            {
                tangents[cellI] = curCf[upperFace];
            }
            else
            {
                label patchID = -upperFace - 1;
                tangents[cellI] = curCf.boundaryField()[patchID][0];
            }

            if (lowerFace > -1)
            {
                tangents[cellI] -= curCf[lowerFace];
            }
            else
            {
                label patchID = -lowerFace - 1;
                tangents[cellI] -= curCf.boundaryField()[patchID][0];
            }

            tangents[cellI] /= mag(tangents[cellI]) + SMALL;
        }

        const volVectorField& C = mesh.C();
        const volVectorField curC(C + refW_ + W_);

        // Correct internal faces mean line position
        const labelList& own = mesh.owner();
        const labelList& nei = mesh.neighbour();
        forAll(curCf, faceI)
        {
            curCf[faceI] =
                HermiteSpline::position
                (
                    curC[own[faceI]],
                    tangents[own[faceI]],
                    curC[nei[faceI]],
                    tangents[nei[faceI]],
                    0
                );
        }

        tCurrentPoints.ref() = this->beamPointData(curCf);

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

        const surfaceVectorField Wf(fvc::interpolate(W_));
        surfaceVectorField curCf(mesh.Cf() + refWf_ + Wf);

        // Calculate cell-centre mean line tangents
        vectorField tangents(mesh.nCells(), vector::zero);
        const labelList& upperNeiFaces = upperNeiCellFaces();
        const labelList& lowerNeiFaces = lowerNeiCellFaces();
        forAll(tangents, cellI)
        {
            label upperFace = upperNeiFaces[cellI];
            label lowerFace = lowerNeiFaces[cellI];

            if (upperFace > -1)
            {
                tangents[cellI] = curCf[upperFace];
            }
            else
            {
                label patchID = -upperFace - 1;

                const labelList& fc =
                    mesh.boundary()[patchID].faceCells();

                forAll(curCf.boundaryField()[patchID], fI)
                {
                    if (fc[fI] == cellI)
                    {
                        tangents[cellI] = curCf.boundaryField()[patchID][fI];
                    }
                }
            }

            if (lowerFace > -1)
            {
                tangents[cellI] -= curCf[lowerFace];
            }
            else
            {
                label patchID = -lowerFace - 1;

                const labelList& fc =
                    mesh.boundary()[patchID].faceCells();

                forAll(curCf.boundaryField()[patchID], fI)
                {
                    if (fc[fI] == cellI)
                    {
                        tangents[cellI] -= curCf.boundaryField()[patchID][fI];
                    }
                }
            }

            tangents[cellI] /= mag(tangents[cellI]) + SMALL;
        }

        const volVectorField& C = mesh.C();
        const volVectorField curC(C + refW_ + W_);

        // Correct internal faces mean line position
        const labelList& own = mesh.owner();
        const labelList& nei = mesh.neighbour();
        forAll(curCf, faceI)
        {
            curCf[faceI] =
                HermiteSpline::position
                (
                    curC[own[faceI]],
                    tangents[own[faceI]],
                    curC[nei[faceI]],
                    tangents[nei[faceI]],
                    0
                );
        }

        tCurrentPoints.ref() = this->beamPointData(curCf, bI);

        return tCurrentPoints;
    }
}


tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamPoints
(
    const label bI
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allPoints(Pstream::nProcs());
        allPoints[Pstream::myProcNo()] = currentBeamPoints(bI);

        // Number of point per processors
        labelList procNPoints(Pstream::nProcs(), 0);
        procNPoints[Pstream::myProcNo()] =
            allPoints[Pstream::myProcNo()].size();
        Pstream::gatherList(procNPoints);
        Pstream::scatterList(procNPoints);

        forAll(allPoints, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                allPoints[procI] =
                    pointField(procNPoints[procI], vector::zero);
            }
        }

        // PstreamBuffers pBufs(Pstream::nonBlocking);


        // Send data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                // if (allPoints[Pstream::myProcNo()].size())
                {
                    // Parallel data exchange
                    OPstream::write
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<const char*>
                        (
                            allPoints[Pstream::myProcNo()].begin()
                        ),
                        allPoints[Pstream::myProcNo()].byteSize()
                    );
                }
            }
        }

        // Pout << "data send" << endl;
        // sleep(5);

        // Receive data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                // if (procNPoints[domainI])
                {
                    // Parallel data exchange
                    IPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<char*>(allPoints[domainI].begin()),
                        allPoints[domainI].byteSize()
                    );
                }
            }
        }

        // Pout << "data received" << endl;
        // sleep(5);

        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (allPoints[domainI].size())
            {
                return tmp<pointField>(new pointField(allPoints[domainI]));
                break;
            }
        }
    }

    return currentBeamPoints(bI);
}


tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentBeamTangents
(
    const label bI
) const
{
    const fvMesh& mesh = this->mesh();
    label nCellZones = mesh.cellZones().size();

    if (nCellZones < 2)
    {
        label nPoints = mesh.nCells() + 1;

        tmp<vectorField> tCurrentTangents
        (
            new vectorField(nPoints, vector::zero)
        );

        volVectorField R0
        (
            IOobject
            (
                "R0",
                runTime().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
             ),
            mesh,
            dimensionedVector("R0", dimLength, vector::zero)
        );
        R0 = mesh.C() + refW_ + W_;


        //R0.boundaryField().evaluateCoupled();
        surfaceVectorField curTangents(fvc::snGrad(R0));
        curTangents /= mag(curTangents);

        if (curTangents.boundaryField()[startPatchIndex()].size())
        {
            curTangents.boundaryFieldRef()[startPatchIndex()] *= -1;
        }
        else
        {
            // Find left processor patch
            forAll(mesh.boundary(), patchI)
            {
                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const labelList& fc =
                        mesh.boundary()[patchI].faceCells();

                    if (fc[0] == 0)
                    {
                        curTangents.boundaryFieldRef()[patchI] *= -1;
                    }
                }
            }
        }

        tCurrentTangents.ref() = this->beamPointData(curTangents);

        return tCurrentTangents;
    }
    else
    {
        const cellZone& cz = mesh.cellZones()[bI];
        label nPoints = cz.size() + 1;

        tmp<vectorField> tCurrentTangents
        (
            new vectorField(nPoints, vector::zero)
        );

        volVectorField R0
        (
            IOobject
            (
                "R0",
                runTime().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("R0", dimLength, vector::zero)
        );
        R0 = mesh.C() + refW_ + W_;
        //R0.boundaryField().evaluateCoupled();
        surfaceVectorField curTangents(fvc::snGrad(R0));
        curTangents /= mag(curTangents);

        // curTangents.boundaryField()[startPatchIndex(bI)] *= -1;

        if (curTangents.boundaryField()[startPatchIndex(bI)].size())
        {
            curTangents.boundaryFieldRef()[startPatchIndex(bI)] *= -1;
        }
        else
        {
            // Find left (start) processor patch
            forAll(mesh.boundary(), patchI)
            {
                if (isA<processorFvPatch>(mesh.boundary()[patchI]))
                {
                    const labelList& fc =
                        mesh.boundary()[patchI].faceCells();

                    forAll(fc, fI)
                    {
                        if (cz.whichCell(fc[fI]) != -1)
                        {
                            if (fc[fI] <= min(cz))
                            {
                                curTangents.boundaryFieldRef()[patchI][fI] *= -1;
                            }
                        }
                    }
                }
            }
        }

        tCurrentTangents.ref() = this->beamPointData(curTangents, bI);

        return tCurrentTangents;
    }
}



tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamTangents
(
    const label bI
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allTangents(Pstream::nProcs());
        allTangents[Pstream::myProcNo()] = currentBeamTangents(bI);

        // Number of point per processors
        labelList procNTangents(Pstream::nProcs(), 0);
        procNTangents[Pstream::myProcNo()] =
            allTangents[Pstream::myProcNo()].size();
        Pstream::gatherList(procNTangents);
        Pstream::scatterList(procNTangents);

        forAll(allTangents, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                allTangents[procI] =
                    pointField(procNTangents[procI], vector::zero);
            }
        }

        // PstreamBuffers pBufs(Pstream::nonBlocking);


        // Send data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                if (allTangents[Pstream::myProcNo()].size())
                {
                    // Parallel data exchange
                    OPstream::write
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<const char*>
                        (
                            allTangents[Pstream::myProcNo()].begin()
                        ),
                        allTangents[Pstream::myProcNo()].byteSize()
                    );
                }
            }
        }

        // Pout << "data send" << endl;
        // sleep(5);

        // Receive data
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (domainI != Pstream::myProcNo())
            {
                if (procNTangents[domainI])
                {
                    // Parallel data exchange
                    IPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        domainI,
                        reinterpret_cast<char*>(allTangents[domainI].begin()),
                        allTangents[domainI].byteSize()
                    );
                }
            }
        }

        // Pout << "data received" << endl;
        // sleep(5);

        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            if (allTangents[domainI].size())
            {
                return tmp<pointField>(new pointField(allTangents[domainI]));
                break;
            }
        }
    }

    return currentBeamTangents(bI);
}


void coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamPointsAndTangents
(
    const label bI,
    vectorField& points,
    vectorField& tangents
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allPoints(Pstream::nProcs());
        allPoints[Pstream::myProcNo()] = currentBeamPoints(bI);

        Pstream::gatherList(allPoints);
        Pstream::scatterList(allPoints);

        List<pointField> allTangents(Pstream::nProcs());
        allTangents[Pstream::myProcNo()] = currentBeamTangents(bI);

        Pstream::gatherList(allTangents);
        Pstream::scatterList(allTangents);

        labelListList glPtIndices(Pstream::nProcs());
        glPtIndices[Pstream::myProcNo()] = globalPointsIndices(bI);

        Pstream::gatherList(glPtIndices);
        Pstream::scatterList(glPtIndices);

        // Set size of points and tangents
        label nPoints = 0;
        forAll(glPtIndices, procI)
        {
            label curMax = max(glPtIndices[procI]);

            if (curMax > nPoints)
            {
                nPoints = curMax;
            }
        }
        nPoints++;

        points.setSize(nPoints);
        tangents.setSize(nPoints);

        // Combine data (points and tangets should be initialised befor)
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            forAll(allPoints[domainI], pI)
            {
                points[glPtIndices[domainI][pI]] = allPoints[domainI][pI];
                tangents[glPtIndices[domainI][pI]] = allTangents[domainI][pI];
            }
        }
    }
    else
    {
        points = currentBeamPoints(bI);
        tangents = currentBeamTangents(bI);
    }
}


void coupledTotalLagNewtonRaphsonBeam::currentGlobalBeamPointsAndTangents
(
    const label bI,
    const labelListList& procSendPoints,
    vectorField& points,
    vectorField& tangents
) const
{
    if (Pstream::parRun())
    {
        List<pointField> allPoints(Pstream::nProcs());
        allPoints[Pstream::myProcNo()] = currentBeamPoints(bI);

        // Pstream::gatherList(allPoints);
        // Pstream::scatterList(allPoints);

        List<pointField> allTangents(Pstream::nProcs());
        allTangents[Pstream::myProcNo()] = currentBeamTangents(bI);

        // Pstream::gatherList(allTangents);
        // Pstream::scatterList(allTangents);

        // labelListList glPtIndices(Pstream::nProcs());
        // glPtIndices[Pstream::myProcNo()] =
        //     localToGlobalBeamPointsAddressing()[bI];
        // glPtIndices[Pstream::myProcNo()] = globalPointsIndices(bI);

        // Pstream::gatherList(glPtIndices);
        // Pstream::scatterList(glPtIndices);



        // if (true){

        // Comunicate number of points which will be send to
        // neighbour processors
        labelListList procNumSendPoints(Pstream::nProcs());

        procNumSendPoints[Pstream::myProcNo()] =
            labelList(Pstream::nProcs(), 0);

        forAll(procNumSendPoints[Pstream::myProcNo()], procI)
        {
            if (procI != Pstream::myProcNo())
            {
                procNumSendPoints[Pstream::myProcNo()][procI] =
                    procSendPoints[procI].size();
            }
            // else
            // {
            //     procNumSendPoints[Pstream::myProcNo()][procI] = 0;
            // }
        }

        Pstream::gatherList(procNumSendPoints);
        Pstream::scatterList(procNumSendPoints);

        // // Number of point per processors
        // labelList procNPoints(Pstream::nProcs(), 0);
        // procNPoints[Pstream::myProcNo()] =
        //     allPoints[Pstream::myProcNo()].size();

        // Pstream::gatherList(procNPoints);
        // Pstream::scatterList(procNPoints);

        forAll(allPoints, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                label curProcNumReceivePoints =
                    procNumSendPoints[procI][Pstream::myProcNo()];
                    // procReceivePoints[procI].size();

                // Pout << curProcNumReceivePoints << endl;

                allPoints[procI] =
                    pointField(curProcNumReceivePoints, vector::zero);
                allTangents[procI] =
                    pointField(curProcNumReceivePoints, vector::zero);

                // allPoints[procI] =
                //     pointField(procNPoints[procI], vector::zero);
                // allTangents[procI] =
                //     pointField(procNPoints[procI], vector::zero);
            }
        }

        // Send and receive point data
        labelListList procReceivePoints_(Pstream::nProcs());
        // procReceivePoints[Pstream::myProcNo()] =
        //     procSendPoints[Pstream::myProcNo()];
        {
            // Send data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (allPoints[Pstream::myProcNo()].size())
                    {
                        // Parallel data exchange
                        OPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<const char*>
                            (
                                procSendPoints[domainI].begin()
                            ),
                            procSendPoints[domainI].byteSize()
                        );
                    }
                }
            }

            // Receive data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    procReceivePoints_[domainI].setSize
                    (
                        procNumSendPoints[domainI][Pstream::myProcNo()]
                    );

                    // if (procNPoints[domainI])
                    {
                        // Parallel data exchange
                        IPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<char*>
                            (
                                procReceivePoints_[domainI].begin()
                            ),
                            procReceivePoints_[domainI].byteSize()
                        );
                    }
                }
            }
        }

        // sleep(5);

        // PstreamBuffers pBufs(Pstream::nonBlocking);

        // Send and receive points data
        {
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (allPoints[Pstream::myProcNo()].size())
                    {
                        vectorField pointsToSend
                        (
                            allPoints[Pstream::myProcNo()],
                            procSendPoints[domainI]
                        );

                        // Parallel data exchange
                        OPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<const char*>
                            (
                                pointsToSend.begin()
                                // allPoints[Pstream::myProcNo()].begin()
                            ),
                            pointsToSend.byteSize()
                            // allPoints[Pstream::myProcNo()].byteSize()
                        );
                    }
                }
            }

            // Receive data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (procNPoints[domainI])
                    {
                        // Parallel data exchange
                        IPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<char*>(allPoints[domainI].begin()),
                            allPoints[domainI].byteSize()
                        );
                    }
                }
            }
        }

        // Send and receive tangents data
        {
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (allPoints[Pstream::myProcNo()].size())
                    {
                        vectorField tangentsToSend
                        (
                            allTangents[Pstream::myProcNo()],
                            procSendPoints[domainI]
                        );

                        OPstream::write
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<const char*>
                            (
                                tangentsToSend.begin()
                                // allTangents[Pstream::myProcNo()].begin()
                            ),
                            tangentsToSend.byteSize()
                            // allTangents[Pstream::myProcNo()].byteSize()
                        );
                    }
                }
            }

            // Receive data
            for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
            {
                if (domainI != Pstream::myProcNo())
                {
                    // if (procNPoints[domainI])
                    {
                        // Parallel data exchange
                        IPstream::read
                        (
                            Pstream::commsTypes::blocking,
                            domainI,
                            reinterpret_cast<char*>
                            (
                                allTangents[domainI].begin()
                            ),
                            allTangents[domainI].byteSize()
                        );
                    }
                }
            }
        }
        // }

        // Combine data
        for (label procI = 0; procI < Pstream::nProcs(); ++procI)
        {
            if (procI != Pstream::myProcNo())
            {
                forAll(allPoints[procI], pI)
                {
                    label locPointIndex = procReceivePoints_[procI][pI];
                    label glPointIndex =
                        localToGlobalBeamPointsAddressing()[bI][procI][locPointIndex];

                    points[glPointIndex] = allPoints[procI][pI];
                    tangents[glPointIndex] = allTangents[procI][pI];

                    // points[glPtIndices[procI][locPointIndex]] = allPoints[procI][pI];
                    // tangents[glPtIndices[procI][locPointIndex]] = allTangents[procI][pI];
                }
            }
            else
            {
                forAll(allPoints[procI], pI)
                {
                    label glPointIndex =
                        localToGlobalBeamPointsAddressing()[bI][procI][pI];

                    points[glPointIndex] = allPoints[procI][pI];
                    tangents[glPointIndex] = allTangents[procI][pI];

                    // points[glPtIndices[procI][pI]] = allPoints[procI][pI];
                    // tangents[glPtIndices[procI][pI]] = allTangents[procI][pI];
                }
            }
        }

        if (false){
        // label nPoints = sum(procNPoints);
        // pointField globalPoints(nPoints, vector::zero);
        // pointField globalTangents(nPoints, vector::zero);
        pointField globalPoints = points;
        pointField globalTangents = tangents;

        // Combine data
        label pointI = 0;
        for (label domainI = 0; domainI < Pstream::nProcs(); ++domainI)
        {
            for (label pI=0; pI<(allPoints[domainI].size()-1); pI++)
            // forAll(allPoints[domainI], pI)
            {
                globalPoints[pointI] = allPoints[domainI][pI];
                globalTangents[pointI] = allTangents[domainI][pI];
                pointI++;
            }

            // Last point
            if (domainI == Pstream::nProcs()-1)
            {
                label pI = allPoints[domainI].size()-1;
                globalPoints[pointI] = allPoints[domainI][pI];
                globalTangents[pointI] = allTangents[domainI][pI];
                pointI++;
            }

            // if (allPoints[domainI].size())
            // {
            //     points = allPoints[domainI];
            //     tangents = allTangents[domainI];
            //     // return tmp<pointField>(new pointField(allPoints[domainI]));
            //     break;
            // }
        }

        globalPoints.resize(pointI);
        globalTangents.resize(pointI);

        points = globalPoints;
        tangents = globalTangents;
        }

        // Pout << "nPoints: " << points.size() << ", " << nPoints << endl;

        // if (false)
        // {
        // // Merge
        // labelList oldToNew;
        // pointField newGlobalPoints;
        // pointField newGlobalTangents;

        // scalar minSegLen = gMin(L().internalField());
        // // scalar minSegLen = gSum(L().internalField());

        // bool hasMerged = mergePoints
        // (
        //     globalPoints,
        //     1e-2*minSegLen, //SMALL,
        //     false,
        //     oldToNew,
        //     newGlobalPoints
        // );

        // if (hasMerged)
        // {
        //     tangents.setSize(newGlobalPoints.size());

        //     // forAll(tangents, pI)
        //     // {
        //     //     label oldIndex = findIndex(oldToNew, pI);
        //     //     tangents[pI] = globalTangents[oldIndex];
        //     // }

        //     forAll(globalTangents, gpI)
        //     {
        //         tangents[oldToNew[gpI]] = globalTangents[gpI];
        //     }

        //     points = newGlobalPoints;
        // }
        // else
        // {
        //     points = globalPoints;
        //     tangents = globalTangents;
        // }
        // }
    }
    else
    {
        points = currentBeamPoints(bI);
        tangents = currentBeamTangents(bI);
    }
}


tmp<vectorField> coupledTotalLagNewtonRaphsonBeam::
currentDisplacementIncrement() const
{
    label nPoints = this->mesh().nCells() + 1;

    tmp<vectorField> tDW
    (
        new vectorField(nPoints, vector::zero)
    );
    vectorField& DW(tDW.ref());

    const surfaceVectorField DWf
    (
        fvc::interpolate(W_)
      - fvc::interpolate(W_.oldTime())
    );

    const vectorField& DWfI = DWf.internalField();

    DW[0] = DWf.boundaryField()[startPatchIndex()][0];
    DW[nPoints-1] = DWf.boundaryField()[endPatchIndex()][0];
    for (label i=0; i<DWfI.size(); i++)
    {
        DW[i+1] = DWfI[i];
    }

    return tDW;
}

tmp<tensorField> coupledTotalLagNewtonRaphsonBeam::
currentRotationIncrement() const
{
    label nPoints = this->mesh().nCells() + 1;

    tmp<tensorField> tDLambda
    (
        new tensorField(nPoints, tensor::zero)
    );
    tensorField& DLambda(tDLambda.ref());

    const surfaceTensorField DLambdaf((Lambda_ & inv(Lambda_.oldTime())));

    const tensorField& DLambdafI(DLambdaf.internalField());

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
