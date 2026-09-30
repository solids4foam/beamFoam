/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Coupling between the beam and a rigidBodyEnd owned by beamFoam.

    rigidBodyCoupling selects the method:
        none         no rigid body (default without a rigidBodyEnd dict)
        partitioned  body and beam solved in turn, with optional coupling
                     iterations (default with a rigidBodyEnd dict)
        monolithic   body translation solved inside the BlockEigen system,
                     in the same Newton iteration as the beam (Phase 1 of
                     codexLogs/monolithicCouplingPlan; rotation not yet)

    The legacy blockEigen* switches are used by the moorFV-driven path and
    cannot be combined with rigidBodyEnd.

\*---------------------------------------------------------------------------*/

#include "coupledTotalLagNewtonRaphsonBeam.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledTotalLagNewtonRaphsonBeam::readRigidBodyEnd()
{
    const dictionary& coeffs = beamProperties();
    const bool hasRigidBodyEnd = coeffs.found("rigidBodyEnd");

    rigidBodyCoupling_ =
        coeffs.getOrDefault<word>
        (
            "rigidBodyCoupling",
            hasRigidBodyEnd ? "partitioned" : "none"
        );

    if
    (
        !List<word>({"none", "partitioned", "monolithic"})
            .found(rigidBodyCoupling_)
    )
    {
        FatalIOErrorInFunction(coeffs)
            << "Invalid rigidBodyCoupling " << rigidBodyCoupling_ << nl
            << "Valid options are none, partitioned and monolithic"
            << exit(FatalIOError);
    }

    const wordList legacySwitches
    ({
        "blockEigenForceCoupling",
        "blockEigenKinematicCoupling"
    });

    bool legacyCouplingOn = false;
    forAll(legacySwitches, i)
    {
        legacyCouplingOn =
            legacyCouplingOn
         || coeffs.getOrDefault<Switch>(legacySwitches[i], false);
    }

    if (!hasRigidBodyEnd)
    {
        if (rigidBodyCoupling_ != "none")
        {
            FatalIOErrorInFunction(coeffs)
                << "rigidBodyCoupling " << rigidBodyCoupling_
                << " needs a rigidBodyEnd dictionary"
                << exit(FatalIOError);
        }

        if (legacyCouplingOn)
        {
            Info<< "Legacy blockEigen coupling switches are set. They are "
                << "used by the moorFV-driven path and do not give a "
                << "monolithic beam / rigid-body solve." << endl;
        }

        return;
    }

    if (rigidBodyCoupling_ == "none")
    {
        FatalIOErrorInFunction(coeffs)
            << "A rigidBodyEnd dictionary is given but rigidBodyCoupling "
            << "is none" << exit(FatalIOError);
    }

    if (legacyCouplingOn)
    {
        FatalIOErrorInFunction(coeffs)
            << "rigidBodyEnd cannot be combined with the legacy "
            << legacySwitches << " switches" << exit(FatalIOError);
    }

    const word defaultPatchName =
        endPatchIndex() >= 0
      ? mesh().boundaryMesh()[endPatchIndex()].name()
      : word("right");

    rigidBodyEndPtr_.reset
    (
        new rigidBodyEnd
        (
            runTime(),
            coeffs.subDict("rigidBodyEnd"),
            defaultPatchName,
            g().value()
        )
    );

    const label patchI =
        mesh().boundaryMesh().findPatchID(rigidBodyEndPtr_().patchName());

    if (patchI < 0)
    {
        FatalErrorInFunction
            << "rigidBodyEnd patch " << rigidBodyEndPtr_().patchName()
            << " not found" << abort(FatalError);
    }

    if (!isA<fixedValueFvPatchVectorField>(W_.boundaryField()[patchI]))
    {
        FatalErrorInFunction
            << "rigidBodyEnd needs a fixedValue W boundary condition on "
            << "patch " << rigidBodyEndPtr_().patchName()
            << abort(FatalError);
    }

    if (W_.boundaryField()[patchI].size() != 1)
    {
        FatalErrorInFunction
            << "rigidBodyEnd patch " << rigidBodyEndPtr_().patchName()
            << " must have exactly one face on this processor"
            << abort(FatalError);
    }

    rigidBodyEndPatchIndex_ = patchI;
    rigidBodyEndPtr_().setAttachmentOffset(W_.boundaryField()[patchI][0]);
    rigidBodyEndPtr_().setAttachmentPoint
    (
        mesh().Cf().boundaryField()[patchI][0] + W_.boundaryField()[patchI][0]
    );

    Info<< "rigidBodyCoupling " << rigidBodyCoupling_ << endl;
}


scalar coupledTotalLagNewtonRaphsonBeam::evolve()
{
    if (rigidBodyEndPtr_.valid())
    {
        if (rigidBodyDataValid_)
        {
            FatalErrorInFunction
                << "rigidBodyEnd cannot be used when rigid-body data is also "
                << "supplied by moorFV" << abort(FatalError);
        }

        if (rigidBodyCoupling_ == "monolithic")
        {
            return evolveMonolithicRigidBodyEnd();
        }

        return evolvePartitionedRigidBodyEnd();
    }

    return evolveBeam();
}


scalar coupledTotalLagNewtonRaphsonBeam::evolvePartitionedRigidBodyEnd()
{
    rigidBodyEnd& body = rigidBodyEndPtr_();

    body.newTimeStep();

    const label patchI = rigidBodyEndPatchIndex_;

    // Length that converts a rotation change to a displacement change
    const scalar armLength = mag(body.arm(vector::zero));

    // Coupling starts from the beam force of the last accepted state
    vector beamForce = body.beamForce();
    vector x = body.displacement();
    vector theta = vector::zero;

    scalar initialResidual = 0;
    label nBeamSolves = 0;

    for (label iter = 1; iter <= body.nCouplingIterations(); ++iter)
    {
        const vector xNewmark = body.newmarkDisplacement(beamForce);
        const vector thetaNewmark = body.newmarkRotation(beamForce, theta);

        // Change relative to the larger of the step and total motion, so the
        // test still works when the body is at rest
        const scalar change =
            (mag(xNewmark - x) + armLength*mag(thetaNewmark - theta))
           /max
            (
                max
                (
                    mag(xNewmark - body.oldDisplacement())
                  + armLength*mag(thetaNewmark),
                    mag(xNewmark)
                ),
                VSMALL
            );

        if (iter > 1)
        {
            Info<< "rigidBodyEnd coupling iteration " << iter - 1
                << ": relative change " << change << endl;

            if (change < body.couplingTolerance())
            {
                break;
            }
        }

        const scalar w = iter > 1 ? body.couplingRelaxation() : 1;
        x = w*xNewmark + (1 - w)*x;
        theta = w*thetaNewmark + (1 - w)*theta;

        W_.boundaryFieldRef()[patchI][0] = rigidBodyEndAttachmentW(x, theta);

        const scalar residual = evolveBeam();

        if (nBeamSolves == 0)
        {
            initialResidual = residual;
        }

        ++nBeamSolves;

        // Force from the beam on the body
        beamForce = -Q_.boundaryField()[patchI][0];
    }

    body.accept(x, theta, beamForce, nBeamSolves);
    body.write();

    Info<< "rigidBodyEnd: displacement " << body.displacement()
        << ", velocity " << body.velocity()
        << ", angular velocity " << body.angularVelocity()
        << ", beam force " << body.beamForce()
        << ", beam solves " << nBeamSolves << endl;

    return initialResidual;
}


scalar coupledTotalLagNewtonRaphsonBeam::evolveMonolithicRigidBodyEnd()
{
    rigidBodyEnd& body = rigidBodyEndPtr_();

    body.newTimeStep();

    // Start from the last accepted state; the Newton iterations of
    // evolveBeam() then update the body and the beam together
    setRigidBodyEndIterate(body.displacement(), vector::zero);

    const scalar residual = evolveBeam();

    body.accept
    (
        rigidBodyEndX_,
        rigidBodyEndTheta_,
        -Q_.boundaryField()[rigidBodyEndPatchIndex_][0],
        iOuterCorr()
    );
    body.write();

    Info<< "rigidBodyEnd: displacement " << body.displacement()
        << ", velocity " << body.velocity()
        << ", angular velocity " << body.angularVelocity()
        << ", beam force " << body.beamForce()
        << ", Newton iterations " << iOuterCorr() << endl;

    return residual;
}


vector coupledTotalLagNewtonRaphsonBeam::rigidBodyEndAttachmentW
(
    const vector& x,
    const vector& theta
) const
{
    const rigidBodyEnd& body = rigidBodyEndPtr_();

    return body.attachmentDisplacement(x, theta) + body.attachmentOffset();
}


void coupledTotalLagNewtonRaphsonBeam::setRigidBodyEndIterate
(
    const vector& x,
    const vector& theta
)
{
    rigidBodyEndX_ = x;
    rigidBodyEndTheta_ = theta;

    W_.boundaryFieldRef()[rigidBodyEndPatchIndex_][0] =
        rigidBodyEndAttachmentW(x, theta);

    // W_ is otherwise unchanged since its last storePrevIter(), so only the
    // attachment value of the previous iterate needs to follow
    W_.storePrevIter();
}


RigidBodyMonolithicCoupling
coupledTotalLagNewtonRaphsonBeam::rigidBodyEndMonolithicCoupling()
{
    const rigidBodyEnd& body = rigidBodyEndPtr_();
    const label patchI = rigidBodyEndPatchIndex_;
    const label faceI = 0;
    const label cellI = mesh().boundary()[patchI].faceCells()[faceI];

    const scalar pDelta =
        1.0/mesh().deltaCoeffs().boundaryField()[patchI][faceI];

    // Attachment face flux per unit boundary displacement: these are the
    // coefficients of the fixedValue boundary increment in the attachment
    // cell's force and moment rows
    const tensor Kb = CQW_.boundaryField()[patchI][faceI]/pDelta;
    const tensor KbM = CMQW_.boundaryField()[patchI][faceI]/pDelta;

    // Boundary contributions on their own give the attachment face force
    // as the beam equations see it, and its dependence on the attachment
    // cell unknowns
    Field<scalarSquareMatrix> dB
    (
        mesh().nCells(), scalarSquareMatrix(6, 0.0)
    );
    Field<scalarSquareMatrix> lB
    (
        mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0)
    );
    Field<scalarSquareMatrix> uB
    (
        mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0)
    );
    Field<scalarRectangularMatrix> sB
    (
        mesh().nCells(), scalarRectangularMatrix(6, 1, 0.0)
    );

    assembleBoundaryConditions(dB, lB, uB, sB);

    const scalarSquareMatrix& dc = dB[cellI];
    const vector Qb(-sB[cellI](0, 0), -sB[cellI](1, 0), -sB[cellI](2, 0));

    const tensor dcW
    (
        dc(0, 0), dc(0, 1), dc(0, 2),
        dc(1, 0), dc(1, 1), dc(1, 2),
        dc(2, 0), dc(2, 1), dc(2, 2)
    );
    const tensor dcTheta
    (
        dc(0, 3), dc(0, 4), dc(0, 5),
        dc(1, 3), dc(1, 4), dc(1, 5),
        dc(2, 3), dc(2, 4), dc(2, 5)
    );

    // Arm from the centre of mass to the attachment point; a change dTheta
    // of the step rotation vector turns the body by T dTheta and moves the
    // attachment point by (T dTheta) x r = -Sr T dTheta
    const vector r = body.arm(rigidBodyEndTheta_);
    const tensor Sr = *r;
    const tensor SrT = (Sr & body.rotationJacobian(rigidBodyEndTheta_));

    RigidBodyMonolithicCoupling coupling;

    coupling.active = true;
    coupling.attachmentCell = cellI;
    coupling.beamWRowCoeff = Kb;
    coupling.beamThetaRowCoeff = KbM;

    // Force balance, with the force from the beam on the body equal to minus
    // the attachment face force:
    //     R = m a + c v - externalForce - m g + Qb
    coupling.bodyCoeff = body.translationalJacobian()*tensor::I + Kb;
    coupling.bodyWCoeff = dcW;
    coupling.bodyThetaCoeff = dcTheta;

    coupling.bodySource = -body.translationalResidual(rigidBodyEndX_, -Qb);

    if (body.rotation())
    {
        coupling.rotationActive = true;

        coupling.beamWRowRotCoeff = -(Kb & SrT);
        coupling.beamThetaRowRotCoeff = -(KbM & SrT);
        coupling.bodyTransRotCoeff = -(Kb & SrT);

        // Moment balance about the centre of mass:
        //     J alpha + omega x (J omega) + cr omega - externalMoment
        //   + r x Qb
        coupling.bodyRotTransCoeff = (Sr & Kb);
        coupling.bodyRotCoeff =
            body.rotationalInertiaJacobian(rigidBodyEndTheta_)
          + ((*Qb) & SrT)
          - (Sr & Kb & SrT);
        coupling.bodyRotWCoeff = (Sr & dcW);
        coupling.bodyRotThetaCoeff = (Sr & dcTheta);

        coupling.bodyRotSource =
            -body.rotationalResidual(rigidBodyEndTheta_, -Qb);
    }

    return coupling;
}


void coupledTotalLagNewtonRaphsonBeam::checkRigidBodyEndJacobian()
{
    const rigidBodyEnd& body = rigidBodyEndPtr_();
    const label patchI = rigidBodyEndPatchIndex_;
    const label cellI = mesh().boundary()[patchI].faceCells()[0];
    const vector x0 = rigidBodyEndX_;
    const vector theta0 = rigidBodyEndTheta_;

    // Strain at the current W; normally updated in updateSolutionVariables()
    auto updateStrain = [&]()
    {
        const surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));
        Gamma_ = (refLambdaf_.T() & ((Lambdaf_.T() & dRdS) - dR0Ds_));
    };

    // Residual of the full system at a body state, with the beam unknowns
    // fixed: beam rows from the complete assembly, body rows from the
    // coupling blocks
    auto evaluate = [&]
    (
        const vector& x,
        const vector& theta,
        vector& beamW,
        vector& beamTheta
    )
    {
        setRigidBodyEndIterate(x, theta);
        updateStrain();
        W_.boundaryFieldRef().updateCoeffs();
        Theta_.boundaryFieldRef().updateCoeffs();
        updateEqnCoefficients();

        Field<scalarSquareMatrix> d
        (
            mesh().nCells(), scalarSquareMatrix(6, 0.0)
        );
        Field<scalarSquareMatrix> l
        (
            mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0)
        );
        Field<scalarSquareMatrix> u
        (
            mesh().nInternalFaces(), scalarSquareMatrix(6, 0.0)
        );
        Field<scalarRectangularMatrix> source
        (
            mesh().nCells(), scalarRectangularMatrix(6, 1, 0.0)
        );

        assembleMatrixCoefficients(d, l, u, source);

        beamW = vector(source[cellI](0, 0), source[cellI](1, 0), source[cellI](2, 0));
        beamTheta = vector(source[cellI](3, 0), source[cellI](4, 0), source[cellI](5, 0));

        return rigidBodyEndMonolithicCoupling();
    };

    vector beamW0, beamTheta0;
    const RigidBodyMonolithicCoupling c0 =
        evaluate(x0, theta0, beamW0, beamTheta0);

    const scalar epsX =
        1e-6*max(mag(x0), 1e-3*mag(mesh().bounds().span()));
    const scalar epsTheta = 1e-7;

    Info<< "rigidBodyEnd Jacobian check (eps " << epsX << " m, " << epsTheta
        << " rad): relative difference between assembled and "
        << "finite-difference columns (entries below 1e-5 of the largest "
        << "body entry in the column count as zero)" << endl;

    const label nColumns = body.rotation() ? 6 : 3;

    for (label col = 0; col < nColumns; ++col)
    {
        const bool rotationColumn = col >= 3;
        const direction i = col % 3;
        const scalar eps = rotationColumn ? epsTheta : epsX;

        vector dx = vector::zero;
        vector dTheta = vector::zero;
        (rotationColumn ? dTheta : dx)[i] = eps;

        vector beamW1, beamTheta1;
        const RigidBodyMonolithicCoupling c1 =
            evaluate(x0 + dx, theta0 + dTheta, beamW1, beamTheta1);

        // Rows are A*dU = source = -residual, so a column of A is
        // -d(source)/d(unknown)
        const vector fdBeamW = -(beamW1 - beamW0)/eps;
        const vector fdBeamTheta = -(beamTheta1 - beamTheta0)/eps;
        const vector fdBody = -(c1.bodySource - c0.bodySource)/eps;
        const vector fdBodyRot = -(c1.bodyRotSource - c0.bodyRotSource)/eps;

        const vector aBeamW
        (
            rotationColumn
          ? c0.beamWRowRotCoeff.col(i)
          : c0.beamWRowCoeff.col(i)
        );
        const vector aBeamTheta
        (
            rotationColumn
          ? c0.beamThetaRowRotCoeff.col(i)
          : c0.beamThetaRowCoeff.col(i)
        );
        const vector aBody
        (
            rotationColumn
          ? c0.bodyTransRotCoeff.col(i)
          : c0.bodyCoeff.col(i)
        );
        const vector aBodyRot
        (
            rotationColumn
          ? c0.bodyRotCoeff.col(i)
          : c0.bodyRotTransCoeff.col(i)
        );

        // Entries much smaller than the body entries count as zero
        const scalar floor = 1e-5*max(mag(aBody), mag(aBodyRot));

        auto relDiff = [floor](const vector& a, const vector& b)
        {
            return mag(a - b)/max(max(mag(a), mag(b)), max(floor, VSMALL));
        };

        Info<< "    column " << col
            << (rotationColumn ? " (rotation)" : " (translation)")
            << ": beam force rows " << relDiff(aBeamW, fdBeamW)
            << " (|col| " << mag(aBeamW) << ")"
            << ", beam moment rows " << relDiff(aBeamTheta, fdBeamTheta)
            << " (|col| " << mag(aBeamTheta) << ")"
            << ", body force rows " << relDiff(aBody, fdBody)
            << " (|col| " << mag(aBody) << ")";

        if (body.rotation())
        {
            Info<< ", body moment rows " << relDiff(aBodyRot, fdBodyRot)
                << " (|col| " << mag(aBodyRot) << ")";
        }

        Info<< endl;
    }

    // Restore the unperturbed state; evolveBeam() recomputes the
    // coefficients before assembling
    setRigidBodyEndIterate(x0, theta0);
    updateStrain();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace beamModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
