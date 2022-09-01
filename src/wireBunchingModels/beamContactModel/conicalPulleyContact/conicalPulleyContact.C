/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

\*---------------------------------------------------------------------------*/

#include "conicalPulleyContact.H"
#include "beamModel.H"
// #include "scalarMatrices.H"
#include "pseudoVector.H"
#include "spinTensor.H"
#include "lagrangeMultipliers.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::conicalPulleyContact::smallFcn_ = 1e-3;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conicalPulleyContact::conicalPulleyContact()
:
    conicalPulleyContactsPtr_(nullptr),
    firstBeam_(-1),
    firstBeamSegment_(-1),
    firstBeamZeta_(0),
    pulleyIndex_(-1),
    pulleyContactCoordinates_(2, vector(0, 0, 0)),
    oldPulleyContactCoordinates_(2, vector(0, 0, 0)),
    delta_(2, GREAT),
    normalDir_(2, vector::zero),
    contactAngle_(2, 0),
    normalContactForce_(2, vector::zero),
    normalContactForceOffset_(2, 0),
    oldNormalContactForce_(2, vector::zero),
    normalContactForceDerivative_(2, tensor::zero),
    normalContactForceDirectionDerivative_(2, tensor::zero),
    normalGap_(2, GREAT),
    oldNormalGap_(2, GREAT),
    frictionalContactMoment_(2, vector::zero),
    frictionalContactForce_(2, vector::zero),
    frictionalContactMomentDerivative_(2, tensor::zero),
    frictionalContactMomentDerivativeOverFcn_(2, tensor::zero),
    frictionalContactForceDerivative_(2, tensor::zero),
    frictionalContactForceDerivativeOverFcn_(2, tensor::zero),
    torsionTheta_(2, 0),
    torsionDTheta_(2, 0),
    oldRM_(2, tensor::I),
    tangentialGap_(2, 0),
    tangentialGapOffset_(2, 0),
    axialTangentialGap_(2, 0),
    axialTangentialGapIncrement_(2, 0),
    axialTangentialGapOffset_(2, 0),
    axialFrictionalContactForce_(2, vector::zero),
    axialFrictionalContactForceDerivative_(2, tensor::zero),
    transversalTangentialGap_(2, 0),
    transversalTangentialGapIncrement_(2, 0),
    transversalTangentialGapOffset_(2, 0),
    transversalFrictionalContactForce_(2, vector::zero),
    transversalFrictionalContactForceDerivativePerDTheta_(2, tensor::zero),
    transversalFrictionalContactForceDerivativePerDW_(2, tensor::zero),
    transversalFrictionalContactMoment_(2, vector::zero),
    transversalFrictionalContactMomentDerivativePerDTheta_(2, tensor::zero),
    transversalFrictionalContactMomentDerivativePerDW_(2, tensor::zero),
    axialBearingFrictionForce_(2, vector::zero),
    normalContactState_(0),
    oldNormalContactState_(0),
    activeFriction_(2, false),
    stick_(2, true),
    unloading_(2, false),
    curTimeIndex_(-1)
{}


Foam::conicalPulleyContact::conicalPulleyContact
(
     const PtrList<conicalPulleyContactListList>& conicalPulleyContacts
)
:
    conicalPulleyContactsPtr_(&conicalPulleyContacts),
    firstBeam_(-1),
    firstBeamSegment_(-1),
    firstBeamZeta_(0),
    pulleyIndex_(-1),
    pulleyContactCoordinates_(2, vector(0, 0, 0)),
    oldPulleyContactCoordinates_(2, vector(0, 0, 0)),
    delta_(2, GREAT),
    normalDir_(2, vector::zero),
    contactAngle_(2, 0),
    normalContactForce_(2, vector::zero),
    normalContactForceOffset_(2, 0),
    oldNormalContactForce_(2, vector::zero),
    normalContactForceDerivative_(2, tensor::zero),
    normalContactForceDirectionDerivative_(2, tensor::zero),
    normalGap_(2, GREAT),
    oldNormalGap_(2, GREAT),
    frictionalContactMoment_(2, vector::zero),
    frictionalContactForce_(2, vector::zero),
    frictionalContactMomentDerivative_(2, tensor::zero),
    frictionalContactMomentDerivativeOverFcn_(2, tensor::zero),
    frictionalContactForceDerivative_(2, tensor::zero),
    frictionalContactForceDerivativeOverFcn_(2, tensor::zero),
    torsionTheta_(2, 0),
    torsionDTheta_(2, 0),
    oldRM_(2, tensor::I),
    tangentialGap_(2, 0),
    tangentialGapOffset_(2, 0),
    axialTangentialGap_(2, 0),
    axialTangentialGapIncrement_(2, 0),
    axialTangentialGapOffset_(2, 0),
    axialFrictionalContactForce_(2, vector::zero),
    axialFrictionalContactForceDerivative_(2, tensor::zero),
    transversalTangentialGap_(2, 0),
    transversalTangentialGapIncrement_(2, 0),
    transversalTangentialGapOffset_(2, 0),
    transversalFrictionalContactForce_(2, vector::zero),
    transversalFrictionalContactForceDerivativePerDTheta_(2, tensor::zero),
    transversalFrictionalContactForceDerivativePerDW_(2, tensor::zero),
    transversalFrictionalContactMoment_(2, vector::zero),
    transversalFrictionalContactMomentDerivativePerDTheta_(2, tensor::zero),
    transversalFrictionalContactMomentDerivativePerDW_(2, tensor::zero),
    axialBearingFrictionForce_(2, vector::zero),
    normalContactState_(0),
    oldNormalContactState_(0),
    activeFriction_(2, false),
    stick_(2, true),
    unloading_(2, false),
    curTimeIndex_(-1)
{}


Foam::conicalPulleyContact::conicalPulleyContact
(
    const PtrList<conicalPulleyContactListList>& conicalPulleyContacts,
    const label fb,
    const label fbSeg,
    const scalar fbZeta,
    const label pulleyIndex,
    const vectorField& pcc,
    const label timeIndex
)
:
    conicalPulleyContactsPtr_(&conicalPulleyContacts),
    firstBeam_(fb),
    firstBeamSegment_(fbSeg),
    firstBeamZeta_(fbZeta),
    pulleyIndex_(pulleyIndex),
    pulleyContactCoordinates_(pcc),
    oldPulleyContactCoordinates_(pcc),
    delta_(2, GREAT),
    normalDir_(2, vector::zero),
    contactAngle_(2, 0),
    normalContactForce_(2, vector::zero),
    normalContactForceOffset_(2, 0),
    oldNormalContactForce_(2, vector::zero),
    normalContactForceDerivative_(2, tensor::zero),
    normalContactForceDirectionDerivative_(2, tensor::zero),
    normalGap_(2, GREAT),
    oldNormalGap_(2, GREAT),
    frictionalContactMoment_(2, vector::zero),
    frictionalContactForce_(2, vector::zero),
    frictionalContactMomentDerivative_(2, tensor::zero),
    frictionalContactMomentDerivativeOverFcn_(2, tensor::zero),
    frictionalContactForceDerivative_(2, tensor::zero),
    frictionalContactForceDerivativeOverFcn_(2, tensor::zero),
    torsionTheta_(2, 0),
    torsionDTheta_(2, 0),
    oldRM_(2, tensor::I),
    // oldTorsionDTheta_(2, 0),
    // oldRM_(2, tensor::zero),
    tangentialGap_(2, 0),
    tangentialGapOffset_(2, 0),
    axialTangentialGap_(2, 0),
    axialTangentialGapIncrement_(2, 0),
    axialTangentialGapOffset_(2, 0),
    axialFrictionalContactForce_(2, vector::zero),
    axialFrictionalContactForceDerivative_(2, tensor::zero),
    transversalTangentialGap_(2, 0),
    transversalTangentialGapIncrement_(2, 0),
    transversalTangentialGapOffset_(2, 0),
    transversalFrictionalContactForce_(2, vector::zero),
    transversalFrictionalContactForceDerivativePerDTheta_(2, tensor::zero),
    transversalFrictionalContactForceDerivativePerDW_(2, tensor::zero),
    transversalFrictionalContactMoment_(2, vector::zero),
    transversalFrictionalContactMomentDerivativePerDTheta_(2, tensor::zero),
    transversalFrictionalContactMomentDerivativePerDW_(2, tensor::zero),
    axialBearingFrictionForce_(2, vector::zero),
    normalContactState_(0),
    oldNormalContactState_(0),
    activeFriction_(2, false),
    stick_(2, true),
    unloading_(2, false),
    curTimeIndex_(timeIndex)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conicalPulleyContact::set
(
    const label fb,
    const label fbSeg,
    const scalar fbZeta,
    const label pulleyIndex,
    const vectorField& pulleyContactCoordinates,
    const label timeIndex
)
{
    if (curTimeIndex_ < timeIndex)
    {
        oldPulleyContactCoordinates_ = pulleyContactCoordinates_;
    }
    
    firstBeam_ = fb;
    firstBeamSegment_ = fbSeg;
    firstBeamZeta_ = fbZeta;

    pulleyIndex_ = pulleyIndex;
    pulleyContactCoordinates_ = pulleyContactCoordinates;

    // Info << conicalPulleyContacts().size() << endl;
}


void Foam::conicalPulleyContact::updateNormalForce
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines,
    const normalContactModel& normalModel,
    const frictionContactModel& frictionModel,
    const scalar contactAngleLimit,
    const bool agumentedLagrangian
)
{
    label bI = firstBeam();
    label segI = firstBeamSegment();
    scalar zeta = firstBeamZeta();

    label nSegments = splines[bI].nSegments();
    
    label pulleyI = pulleyIndex();
    
    bool frictionless = 
        beam.conicalPulleys()[pulleyI].frictionless();
    bool frictionlessFullContact = 
        beam.conicalPulleys()[pulleyI].frictionlessFullContact();

    if (curTimeIndex_ < beam.runTime().timeIndex())
    {
        // this will be don later
        // curTimeIndex_ = beam.runTime().timeIndex();

        vectorField oldOldNormalContactForce =
            oldNormalContactForce_;

        oldNormalContactForce_ = normalContactForce_;
        oldNormalGap_ = normalGap_;
    }

    const vectorField& cylCoord = pulleyContactCoordinates();

    scalar dist = GREAT;
    vector n = vector::one;

    if (pulleyI != -1)
    {
        tensor curRM =
            beam.beamSegmentData
            (
                beam.solutionRMC(),
                bI,
                segI
            )()[0];

	for (label cpI=0; cpI<2; cpI++)
	{
            vector curContactPoint = 
                splines[bI].position(segI, zeta);

	    vector curNeiContactPoint = 
 	        beam.conicalPulleys()[pulleyI].position(cylCoord[cpI]);

	    dist = mag(curContactPoint - curNeiContactPoint);

	    n = curContactPoint - curNeiContactPoint;
	    n /= mag(n) + SMALL;

	    scalar gap = dist - beam.R(bI);

	    normalGap_[cpI] = gap;

	    delta()[cpI] = cylCoord[cpI].x(); //*::sqrt(2.0);
	    // delta()[cpI] = dist;

            // vector nCorr = n - T[segI]*(T[segI] & n);
            // nCorr /= mag(nCorr) + SMALL;

            // Subtract axial component
            vector N = (curRM.T() & n);
            N.x() = 0;
            // Info << N << endl;
            n = (curRM & N);
	    n /= mag(n) + SMALL;

            normalDir_[cpI] = n;

            // if (segI == 985 && pulleyI == 0)
            // {
            //     Info << beam.iOuterCorr() << ", " << gap << endl;
            // }

            if (isA<lagrangeMultipliers>(normalModel))
            {
                // if
                // (
                //     (oldNormalContactForce()[cpI]) < SMALL
                //  && (normalContactForce()[cpI]) < SMALL
                // )
                // {
                //     normalContactForce()[cpI] =
                //         normalModel.contactForce
                //         (
                //             gap,
                //             normalContactForceOffset()[cpI]
                //         )*n;
                // }
                // else
                // {
                //     scalar curNormalContactForce =
                //         sign(normalContactForce()[cpI]&n)
                //        *mag(normalContactForce()[cpI]);
                    
                //     normalContactForce()[cpI] = curNormalContactForce*n;
                // }

                if ((normalContactForce()[cpI] & n) > SMALL)
                {
                    // Update normal contact force direction
                     scalar curNormalContactForce =
                        sign(normalContactForce()[cpI] & n)
                       *mag(normalContactForce()[cpI]);

                    normalContactForce()[cpI] = curNormalContactForce*n;
                }
                else
                {
                    // if (mag(normalContactForce()[cpI]) > SMALL)
                    // {
                    //     Info << segI << ", "
                    //          << (normalContactForce()[cpI] & n) << endl;
                    // }
                    
                    // Normal contact force is zero or positive
                    normalContactForce()[cpI] = vector::zero;
                }

                // scalar curNormalContactForce =
                // (
                //     normalContactForce()[cpI]
                //   & normalDirection()[cpI]
                // );
                    
                // if (gap < 0)
                // {
                //     Info << segI << ", " << beam.iOuterCorr() << ", "
                //          << gap << ", " << curNormalContactForce << ", "
                //          << normalDirection()[cpI] << endl;
                // }
                // else if (curNormalContactForce > SMALL)
                // {
                //     Info << "_" << segI << ", " << beam.iOuterCorr() << ", "
                //          << gap << ", " << curNormalContactForce << ", "
                //          << normalDirection()[cpI] << endl;
                // }
            }
            else
            {
                normalContactForce()[cpI] =
                    normalModel.contactForce
                    (
                        gap,
                        normalContactForceOffset()[cpI]
                    )*n;

                normalContactForceDerivative()[cpI] =
                    normalModel.contactForceDerivative
                    (
                        gap,
                        normalContactForceOffset()[cpI]
                    )*sqr(n);

                if (agumentedLagrangian)
                {
                    // if ((gap < -1e-4*beam.R().value()) || (gap > 0))
                    if
                    (
                        (beam.iOuterCorr() > 2)
                     && (mag(normalContactForce()[cpI]) > SMALL)
                    )
                    {
                        scalar corr =
                        (    
                            mag(normalContactForce()[cpI])
                          - normalContactForceOffset()[cpI]
                        );

                        normalContactForceOffset()[cpI] += corr;

                        normalContactForce()[cpI] =
                            normalModel.contactForce
                            (
                                gap,
                                normalContactForceOffset()[cpI]
                            )*n;

                        normalContactForceDerivative()[cpI] =
                            normalModel.contactForceDerivative
                            (
                                gap,
                                normalContactForceOffset()[cpI]
                            )*sqr(n);
                    }
                    else if (mag(normalContactForceOffset()[cpI]) > SMALL)
                    {
                        // normalContactForceOffset()[cpI] = 0;
                    }
                }
            }

            vector t = curNeiContactPoint;

            if (cpI == 0)
            {
                t -= beam.conicalPulleys()[pulleyI].negConeApex();   
            }
            else
            {
                t -= beam.conicalPulleys()[pulleyI].posConeApex();
            }
            t /= mag(t) + SMALL;

            normalContactForceDirectionDerivative()[cpI] =
                mag(normalContactForce()[cpI])
               *((I-sqr(t))&(I-sqr(n)))/delta()[cpI];
	}
    }
    else
    {
        normalContactForce() = vector::zero;
        normalContactForceDerivative() = tensor::zero;
        normalContactForceDirectionDerivative() = tensor::zero;
    }
}


void Foam::conicalPulleyContact::updateTangentialGap
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines,
    const normalContactModel& normalModel,
    const frictionContactModel& frictionModel,
    const scalar contactAngleLimit,
    const bool agumentedLagrangian
)
{
    label bI = firstBeam();
    label segI = firstBeamSegment();
    scalar zeta = firstBeamZeta();

    label nSegments = splines[bI].nSegments();
    
    label pulleyI = pulleyIndex();
    
    bool frictionless = 
        beam.conicalPulleys()[pulleyI].frictionless();
    bool frictionlessFullContact = 
        beam.conicalPulleys()[pulleyI].frictionlessFullContact();

    if (curTimeIndex_ < beam.runTime().timeIndex())
    {
        curTimeIndex_ = beam.runTime().timeIndex();

        // vectorField oldOldNormalContactForce =
        //     oldNormalContactForce_;
        // oldNormalContactForce_ = normalContactForce_;
        // oldNormalGap_ = normalGap_;

        oldNormalContactState_ = normalContactState_;
        normalContactState_ = 0;
	for (label cpI=0; cpI<2; cpI++)
	{
            if (mag(oldNormalContactForce_[cpI]) > smallFcn_)
            {
                normalContactState_++;
            }
        }

        label normalContactStateDiff =
            normalContactState_
          - oldNormalContactState_;

        // if (normalContactStateDiff == 1) // normal contact state changed
        if (false)
        {
            // activeFriction_ = false;

            axialTangentialGap_ = 0;
            axialTangentialGapIncrement_ = 0;
            axialTangentialGapOffset_ = 0;

            // axialFrictionalContactForce_ = vector::zero;
            // axialFrictionalContactForceDerivative_ = tensor::zero;


            
            // reset tangential gaps

            transversalTangentialGap_ = 0;
            transversalTangentialGapIncrement_ = 0;
            transversalTangentialGapOffset_ = 0;

            // transversalFrictionalContactForce_ = vector::zero;
            // transversalFrictionalContactForceDerivativePerDTheta_ =
            //     tensor::zero;
            // transversalFrictionalContactForceDerivativePerDW_ =
            //     tensor::zero;

            // transversalFrictionalContactMoment_ = vector::zero;
            // transversalFrictionalContactMomentDerivativePerDTheta_ =
            //     tensor::zero;
            // transversalFrictionalContactMomentDerivativePerDW_ =
            //     tensor::zero;
        }

        // Activate/deactivate friction
	for (label cpI=0; cpI<2; cpI++)
	{
            scalar maxFrictionForce =
                frictionModel.maxFrictionForce
                (
                    mag(oldNormalContactForce_[cpI])
                );

            if (!activeFriction_[cpI])
            {
                if
                (
                    (mag(oldNormalContactForce_[cpI]) > smallFcn_)
                 && (maxFrictionForce > SMALL)
                 && !frictionless
                )
                {
                    activeFriction_[cpI] = true;
                    stick_[cpI] = true;

                    // Contact inception: set initial RM
                    oldRM_[cpI] =
                        beam.beamSegmentData
                        (
                            beam.solutionRMC(),
                            // beam.solutionRMC().oldTime(),
                            bI,
                            segI
                        )()[0];

                    torsionTheta_[cpI] = 0;
                    torsionDTheta_[cpI] = 0;

                    // tangentialGap_[cpI] = 0;
                    // tangentialGapOffset_[cpI] = 0;

                    axialTangentialGap_[cpI] = 0;
                    axialTangentialGapIncrement_[cpI] = 0;
                    axialTangentialGapOffset_[cpI] = 0;

                    transversalTangentialGap_[cpI] = 0;
                    transversalTangentialGapIncrement_[cpI] = 0;
                    transversalTangentialGapOffset_[cpI] = 0;
                }
            }
            else
            {
                if( mag(oldNormalContactForce_[cpI]) < smallFcn_ )
                {
                    activeFriction_[cpI] = false;
                    stick_[cpI] = false;
		    
                    oldRM_[cpI] =
                        beam.beamSegmentData
                        (
                            beam.solutionRMC(),
                            // beam.solutionRMC().oldTime(),
                            bI,
                            segI
                        )()[0];

                    torsionTheta_[cpI] = 0;
                    torsionDTheta_[cpI] = 0;

                    tangentialGap_[cpI] = 0;
                    tangentialGapOffset_[cpI] = 0;

                    axialTangentialGap_[cpI] = 0;
                    axialTangentialGapIncrement_[cpI] = 0;
                    axialTangentialGapOffset_[cpI] = 0;

                    transversalTangentialGap_[cpI] = 0;
                    transversalTangentialGapIncrement_[cpI] = 0;
                    transversalTangentialGapOffset_[cpI] = 0;

                    axialFrictionalContactForce_ = vector::zero;
                    axialFrictionalContactForceDerivative_ = tensor::zero;

                    transversalFrictionalContactForce_ = vector::zero;
                    transversalFrictionalContactForceDerivativePerDTheta_ =
                        tensor::zero;
                    transversalFrictionalContactForceDerivativePerDW_ =
                        tensor::zero;

                    transversalFrictionalContactMoment_ = vector::zero;
                    transversalFrictionalContactMomentDerivativePerDTheta_ =
                        tensor::zero;
                    transversalFrictionalContactMomentDerivativePerDW_ =
                        tensor::zero;
                }
                else
                {
                    if (mag(beam.U(bI)) < SMALL)
                    {
                        axialTangentialGap_[cpI] +=
                            axialTangentialGapIncrement_[cpI];

                        // axialTangentialGapIncrement_[cpI] = 0;

                        transversalTangentialGap_[cpI] +=
                            transversalTangentialGapIncrement_[cpI];
                        torsionTheta_[cpI] += torsionDTheta_[cpI];

                        // transversalTangentialGapIncrement_[cpI] = 0;
                        
                        oldRM_[cpI] =
                            beam.beamSegmentData
                            (
                                beam.solutionRMC(),
                                // beam.solutionRMC().oldTime(),
                                bI,
                                segI
                            )()[0];
                        
                        // Axial tangential gap
                        axialTangentialGapOffset_[cpI] =
                            frictionModel.newGapOffset
                            (
                                axialTangentialGap_[cpI],
                                axialTangentialGapIncrement_[cpI],
                                mag(oldNormalContactForce_[cpI]),
                                axialTangentialGapOffset_[cpI]
                            );
                        axialTangentialGapIncrement_[cpI] = 0;

                        // Transversal tangential gap
                        transversalTangentialGapOffset_[cpI] =
                            frictionModel.newGapOffset
                            (
                                transversalTangentialGap_[cpI],
                                transversalTangentialGapIncrement_[cpI],
                                mag(oldNormalContactForce_[cpI]),
                                transversalTangentialGapOffset_[cpI],
                                false
                            );
                        transversalTangentialGapIncrement_[cpI] = 0;
                    }
                    else
                    {
                        // activeFriction_[cpI] = true;
                        // stick_[cpI] = true;

                        // Be carefull with this term for steady-state
                        // (Eulerian) simulations
                        // // Contact inception: set initial RM
                        // oldRM_[cpI] =
                        //     beam.beamSegmentData
                        //     (
                        //         beam.solutionRMC(),
                        //         // beam.solutionRMC().oldTime(),
                        //         bI,
                        //         segI
                        //     )()[0];
                        
                        torsionTheta_[cpI] = 0;
                        torsionDTheta_[cpI] = 0;

                        // tangentialGap_[cpI] = 0;
                        // tangentialGapOffset_[cpI] = 0;

                        axialTangentialGap_[cpI] = 0;
                        axialTangentialGapIncrement_[cpI] = 0;
                        axialTangentialGapOffset_[cpI] = 0;

                        transversalTangentialGap_[cpI] = 0;
                        transversalTangentialGapIncrement_[cpI] = 0;
                        transversalTangentialGapOffset_[cpI] = 0;
                    }
                }
            }
        }
    }

    // const vectorField& cylCoord = pulleyContactCoordinates();



    // Tangential contact
    if (pulleyI != -1 && !frictionless)
    {
        tensor curRM =
            beam.beamSegmentData
            (
                beam.solutionRMC(),
                bI,
                segI
            )()[0];

        vector DW =
            beam.beamSegmentData
            (
                beam.solutionDW(),
                bI,
                segI
            )()[0];

        // vector W =
        //     beam.beamSegmentData
        //     (
        //         beam.solutionW(),
        //         bI,
        //         segI
        //     )()[0];

        // vector oldW =
        //     beam.beamSegmentData
        //     (
        //         beam.solutionW().oldTime(),
        //         bI,
        //         segI
        //     )()[0];

        // vector DW = W - oldW;

        // vector T =
        //     beam.beamSegmentData(beam.solutionCellTangent(), bI, segI)()[0];

        // vector T = splines[bI].dRdS()[segI];
        // T /= mag(T) + SMALL;

        vector T(1, 0, 0); // reference configuration
        T = (curRM & T); // current configuration
        T /= mag(T) + SMALL;

        // vector TT = (curRM.T() & T);
        // TT.y() = 0;
        // TT.z() = 0;
        // T = (curRM & TT);
        // T /= mag(T) + SMALL;

        // vector Tm = (curRM.T() & T);
        // Tm = (oldRM_[0] & Tm);
        // Tm += T;
        // Tm /= mag(Tm) + SMALL;

        // T0 = T;
        // vector T0(1, 0, 0);
        // T0 = (oldRM_[0] & T0);
        // T0 /= mag(T0) + SMALL;

        for (label cpI=0; cpI<2; cpI++)
        {
            if
            (
                activeFriction_[cpI]
             // && (mag(normalContactForce_[cpI]) > smallFcn_)
            )
            {
                if (true)
                // if (beam.iOuterCorr() > 0)
                {
                    vector tangentialGapIncrement = vector::zero;

                    if (mag(beam.U(bI)) > SMALL)
                    {
                        tangentialGapIncrement =
                           -beam.U(bI)*beam.runTime().deltaT().value()*T;

                        // Info << beam.U(bI) << endl;
                    }
                    else
                    {
                        tangentialGapIncrement = -DW;
                    }

                    // note - sign

                    // tangentialGapIncrement -=
                    //     normalDir_[cpI]*(normalDir_[cpI] & DW);

                    vector cylContactPointCoords =
                        pulleyContactCoordinates()[cpI];

                    // scalar rMin =
                    //     beam.conicalPulleys()[pulleyI]
                    //    .fullContactPointRadius(beam.R(bI));

                    // if (cylContactPointCoords.x()<rMin)
                    // {
                    //     cylContactPointCoords.x() = rMin;
                    // }

                    vector locPosition = 
                        beam.conicalPulleys()[pulleyI].position
                        (
                            cylContactPointCoords
                            // pulleyContactCoordinates()[cpI]
                        )
                      - beam.conicalPulleys()[pulleyI].origin();

                    scalar Omega = //0;
                        beam.conicalPulleys()[pulleyI]
                       .angularVelocity();
                       //  beam.conicalPulleys()[pulleyI]
                       // .angularVelocity(beam.R(bI));

                    // Info << "Omega: " << Omega << endl;

                    vector tangentialGapIncrementPulleyRot = // note + sign
                        (beam.conicalPulleys()[pulleyI].axis() ^ locPosition)
                        *Omega*beam.runTime().deltaT().value();
                    
                    // vector tangentialGapIncrement =
                    //     beam.conicalPulleys()[pulleyI].position
                    //     (
                    //         pulleyContactCoordinates()[cpI]
                    //     )
                    //   - beam.conicalPulleys()[pulleyI].position
                    //     (
                    //         oldPulleyContactCoordinates_[cpI]
                    //       // + vector(0, DPhi, 0)
                    //     );

                    if
                    (
                        // true
                           (mag(oldNormalContactForce_[0]) > smallFcn_)
                        && (mag(oldNormalContactForce_[1]) > smallFcn_)
                    )
                    {
                        axialTangentialGap_[cpI] = 0;
                        axialTangentialGapIncrement_[cpI] = 0;
                        // axialFrictionalContactForce_[cpI] = vector::zero;
                        // axialFrictionalContactForceDerivative_[cpI] =
                        //     tensor::zero;
                    }
                    else
                    {
                        // Axial friction force
                        scalar newAxilaTangentialGapIncrement =
                            (T  & tangentialGapIncrement)
                          + (T  & tangentialGapIncrementPulleyRot);

                        // Info << setI << ", "
                        //      << (T & tangentialGapIncrement) << ", "
                        //      << (T & tangentialGapIncrementPulleyRot) << ", "
                        //      << newAxilaTangentialGapIncrement << endl;

                        // axialTangentialGapIncrement_[cpI] =
                        //     (T & tangentialGapIncrement)
                        //   + (T & tangentialGapIncrementPulleyRot);

                        axialTangentialGapIncrement_[cpI] +=
                            // frictionModel.axialRelaxationFactor()
                            (
                                newAxilaTangentialGapIncrement
                              - axialTangentialGapIncrement_[cpI]
                            );

                        // Info << segI << ", "
                        //      << (T & tangentialGapIncrement) << ", "
                        //      << (T & tangentialGapIncrementPulleyRot) << ", "
                        //      << axialTangentialGapIncrement_[cpI] << endl;
                    
                        // Add increment due to wire convection
                        // if (mag(beam.U(bI)) > SMALL)
                        if (false)
                        {
                            label neiSegI = segI;

                            if(beam.U(bI) > SMALL)
                            {
                                if (segI>0)
                                {
                                    neiSegI--;
                                }
                            }
                            else
                            {
                                if (segI<(nSegments-1))
                                {
                                    neiSegI++;
                                }
                            }
                            
                            if (neiSegI!=segI)
                            {
                                const conicalPulleyContact& neiSegContact =
                                    conicalPulleyContacts()[bI][neiSegI][pulleyI];

                                scalar sPN = 
                                    splines[bI].segLength(segI)/2
                                    + splines[bI].segLength(neiSegI)/2;
                                
                                scalar cd =
                                    mag(beam.U(bI))
                                    *beam.runTime().deltaT().value();

                                scalar neiGap =
                                    neiSegContact.axialTangentialGap()[cpI];
                                scalar ownGap =
                                    axialTangentialGap_[cpI];
                                scalar gapCorr = neiGap-ownGap;

                                if (cd<sPN)
                                {
                                    scalar upwindGap =
                                        neiGap
                                      + ((ownGap - neiGap)/(sPN))*(sPN-cd);

                                    gapCorr = upwindGap - ownGap;
                                }
                                else
                                {
                                    Info << "xxxxxxxxxxxxx"
                                        << endl;
                            }

                                axialTangentialGapIncrement_[cpI] += gapCorr;
                            }
                        }

                        // scalar curAxialTangentialGap =
                        //     axialTangentialGap_[cpI]
                        //   + axialTangentialGapIncrement_[cpI];

                        // scalar newAxialFrictionalContactForce =
                        //     frictionModel.contactForce
                        //     (
                        //         curAxialTangentialGap,
                        //         mag(normalContactForce_[cpI]) + smallFcn_,
                        //         axialTangentialGapOffset_[cpI]
                        //     );

                        // axialFrictionalContactForce_[cpI] +=
                        //     frictionModel.axialRelaxationFactor()
                        //    *(
                        //        newAxialFrictionalContactForce
                        //      - (axialFrictionalContactForce_[cpI] & T)
                        //     )*T;

                        // scalar curMaxFrictionForce =
                        //     frictionModel.maxFrictionForce
                        //     (
                        //         mag(normalContactForce_[cpI])
                        //     );
                        
                        // Info << segI << ", " << curAxialTangentialGap << ", "
                        //      << mag(normalContactForce_[cpI]) << ", "
                        //      << (axialFrictionalContactForce_[cpI] & T) << ", "
                        //      << curMaxFrictionForce << ", "
                        //      << axialTangentialGapOffset_[cpI] << endl;

                        // axialFrictionalContactForceDerivative_[cpI] =
                        //    -frictionModel.contactForceDerivative
                        //     (
                        //         curAxialTangentialGap,
                        //         mag(normalContactForce_[cpI]) + smallFcn_,
                        //         axialTangentialGapOffset_[cpI]
                        //     )*sqr(T);
                        
                        // axialFrictionalContactForceDerivative_[cpI] +=
                        //     frictionModel.contactForceDerivative
                        //     (
                        //         curAxialTangentialGap,
                        //         mag(normalContactForce_[cpI]) + smallFcn_,
                        //         axialTangentialGapOffset_[cpI]
                        //     )
                        //    *Omega*beam.runTime().deltaT().value()
                        //    *(
                        //        (
                        //            sqr(T)
                        //            & spinTensor
                        //            (
                        //                beam.conicalPulleys()[pulleyI].axis()
                        //                )
                        //            )
                        //        & (I-sqr(normalDir_[cpI]))
                        //        );
                        
                        // axialFrictionalContactForceDerivative_[cpI] +=
                        //     frictionModel.contactForceDerivativeOverFcn
                        //     (
                        //         curAxialTangentialGap,
                        //         mag(normalContactForce_[cpI]) + smallFcn_,
                        //         axialTangentialGapOffset_[cpI],
                        //         false
                        //     )
                        //    *normalModel.contactForceDerivative
                        //     (
                        //         normalGap_[cpI],
                        //         0 //normalContactForceOffset()
                        //     )
                        //    *(T*normalDir_[cpI]);

                        // Info << "-x2" << endl;
                    }

                    // Transversal fricion force

                    vector R0 = -normalDir_[cpI];
                    R0 -= T*(T&R0);
                    R0 /= mag(R0) + SMALL;
                    // transversal friction force postive direction
                    vector u0 = (T ^ R0);

                    // scalar curTransTangentialGapIncrement = 0;
                    //   (u0 & ((I-sqr(T)) & tangentialGapIncrement))
                    // + (u0 & ((I-sqr(T)) & tangentialGapIncrementPulleyRot));

                    // transversalTangentialGapIncrement_[cpI] = //0;
                    //     curTransTangentialGapIncrement;

                    tensor relRM = (curRM & oldRM_[cpI].T());
                    vector relTheta = pseudoVector(relRM);
                    torsionDTheta_[cpI] =
                        -(oldRM_[cpI].T() & relTheta).x(); // note - sign

                    // torsionDTheta_[cpI] = (T & relTheta); // note - sign
                    // curTransTangentialGapIncrement +=
                    //     beam.R(bI)*torsionDTheta_[cpI];
                    
                    scalar newTransversalTangentialGapIncrement =
                        beam.R(bI)*torsionDTheta_[cpI];

                    if
                    (
                        (mag(oldNormalContactForce_[0]) > smallFcn_)
                     && (mag(oldNormalContactForce_[1]) > smallFcn_)
                    )
                    {
                        // Do nothing
                    }
                    else
                    {
                        newTransversalTangentialGapIncrement +=
                            (u0 & ((I-sqr(T)) & tangentialGapIncrement))
                          + (u0 & ((I-sqr(T)) & tangentialGapIncrementPulleyRot));
                    }

                    transversalTangentialGapIncrement_[cpI] +=
                        // frictionModel.transversalRelaxationFactor()
                        (
                            newTransversalTangentialGapIncrement
                          - transversalTangentialGapIncrement_[cpI]
                        );

                    // Add increment due to wire convection
                    // if (mag(beam.U(bI)) > SMALL)
                    if (false)
                    {
                        label neiSegI = segI;

                        if(beam.U(bI) > SMALL)
                        {
                            if (segI>0)
                            {
                                neiSegI--;
                            }
                        }
                        else
                        {
                            if (segI<(nSegments-1))
                            {
                                neiSegI++;
                            }
                        }

                        if (neiSegI!=segI)
                        {
                            const conicalPulleyContact& neiSegContact =
                                conicalPulleyContacts()[bI][neiSegI][pulleyI];

                            scalar sPN = 
                                splines[bI].segLength(segI)/2
                              + splines[bI].segLength(neiSegI)/2;

                            scalar cd =
                                mag(beam.U(bI))
                               *beam.runTime().deltaT().value();

                            scalar neiGap =
                                neiSegContact.transversalTangentialGap()[cpI];
                            scalar ownGap =
                                transversalTangentialGap_[cpI];
                            scalar gapCorr = neiGap-ownGap;

                            if (cd<sPN)
                            {
                                scalar upwindGap =
                                    neiGap
                                    + ((ownGap - neiGap)/(sPN))*(sPN-cd);
                                
                                gapCorr = upwindGap - ownGap;
                            }
                            else
                            {
                                Info << "xxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
                            }
                            
                            transversalTangentialGapIncrement_[cpI] += gapCorr;
                        }
                    }

                    // scalar curTransversalTangentialGap =
                    //     transversalTangentialGap_[cpI]
                    //   + transversalTangentialGapIncrement_[cpI];
                }
                else
                {
                    // Info << "Predictor: " << beam.iOuterCorr() << endl;
                }
            }
            else
            {
                // frictionalContactMoment_[cpI] = vector::zero;
                // frictionalContactMomentDerivative_[cpI] = tensor::zero;
                // frictionalContactMomentDerivativeOverFcn_[cpI] = tensor::zero;
                // frictionalContactForce_[cpI] = vector::zero;
                // frictionalContactForceDerivative_[cpI] = tensor::zero;
                // frictionalContactForceDerivativeOverFcn_[cpI] = tensor::zero;

                axialFrictionalContactForce_[cpI] = vector::zero;
                axialFrictionalContactForceDerivative_[cpI] = tensor::zero;

                transversalFrictionalContactForce_[cpI] = vector::zero;
                transversalFrictionalContactForceDerivativePerDTheta_[cpI] =
                    tensor::zero;
                transversalFrictionalContactForceDerivativePerDW_[cpI] =
                    tensor::zero;
                
                transversalFrictionalContactMoment_[cpI] = vector::zero;
                transversalFrictionalContactMomentDerivativePerDTheta_[cpI] =
                    tensor::zero;
                transversalFrictionalContactMomentDerivativePerDW_[cpI] =
                    tensor::zero;
            }
        }
    }
}


void Foam::conicalPulleyContact::updateFrictionForceAndMoment
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines,
    const normalContactModel& normalModel,
    const frictionContactModel& frictionModel,
    const scalar contactAngleLimit,
    const bool agumentedLagrangian
)
{
    label bI = firstBeam();
    label segI = firstBeamSegment();
    scalar zeta = firstBeamZeta();

    label nSegments = splines[bI].nSegments();
    
    label pulleyI = pulleyIndex();
    
    bool frictionless = 
        beam.conicalPulleys()[pulleyI].frictionless();
    
    bool frictionlessFullContact = 
        beam.conicalPulleys()[pulleyI].frictionlessFullContact();

    // if (curTimeIndex_ < beam.runTime().timeIndex())
    // {
    //     curTimeIndex_ = beam.runTime().timeIndex();

    //     // vectorField oldOldNormalContactForce =
    //     //     oldNormalContactForce_;
    //     // oldNormalContactForce_ = normalContactForce_;
    //     // oldNormalGap_ = normalGap_;

    //     oldNormalContactState_ = normalContactState_;
    //     normalContactState_ = 0;
    //     for (label cpI=0; cpI<2; cpI++)
    //     {
    //         if (mag(oldNormalContactForce_[cpI]) > smallFcn_)
    //         {
    //             normalContactState_++;
    //         }
    //     }

    //     label normalContactStateDiff =
    //         normalContactState_
    //       - oldNormalContactState_;

    //     // if (normalContactStateDiff == 1) // normal contact state changed
    //     if (false)
    //     {
    //         // activeFriction_ = false;

    //         axialTangentialGap_ = 0;
    //         axialTangentialGapIncrement_ = 0;
    //         axialTangentialGapOffset_ = 0;

    //         axialFrictionalContactForce_ = vector::zero;
    //         axialFrictionalContactForceDerivative_ = tensor::zero;

    //         // reset tangential gaps

    //         transversalTangentialGap_ = 0;
    //         transversalTangentialGapIncrement_ = 0;
    //         transversalTangentialGapOffset_ = 0;

    //         transversalFrictionalContactForce_ = vector::zero;
    //         transversalFrictionalContactForceDerivativePerDTheta_ =
    //             tensor::zero;
    //         transversalFrictionalContactForceDerivativePerDW_ =
    //             tensor::zero;

    //         transversalFrictionalContactMoment_ = vector::zero;
    //         transversalFrictionalContactMomentDerivativePerDTheta_ =
    //             tensor::zero;
    //         transversalFrictionalContactMomentDerivativePerDW_ =
    //             tensor::zero;
    //     }

    //     // Activate/deactivate friction
    //     for (label cpI=0; cpI<2; cpI++)
    //     {
    //         scalar maxFrictionForce =
    //             frictionModel.maxFrictionForce
    //             (
    //                 mag(oldNormalContactForce_[cpI])
    //             );

    //         if (!activeFriction_[cpI])
    //         {
    //             if
    //             (
    //                 (mag(oldNormalContactForce_[cpI]) > smallFcn_)
    //              && (maxFrictionForce > SMALL)
    //              && !frictionless
    //             )
    //             {
    //                 activeFriction_[cpI] = true;
    //                 stick_[cpI] = true;

    //                 // Contact inception: set initial RM
    //                 oldRM_[cpI] =
    //                     beam.beamSegmentData
    //                     (
    //                         beam.solutionRMC(),
    //                         // beam.solutionRMC().oldTime(),
    //                         bI,
    //                         segI
    //                     )()[0];

    //                 torsionTheta_[cpI] = 0;
    //                 torsionDTheta_[cpI] = 0;

    //                 // tangentialGap_[cpI] = 0;
    //                 // tangentialGapOffset_[cpI] = 0;

    //                 axialTangentialGap_[cpI] = 0;
    //                 axialTangentialGapIncrement_[cpI] = 0;
    //                 axialTangentialGapOffset_[cpI] = 0;

    //                 transversalTangentialGap_[cpI] = 0;
    //                 transversalTangentialGapIncrement_[cpI] = 0;
    //                 transversalTangentialGapOffset_[cpI] = 0;
    //             }
    //         }
    //         else
    //         {
    //             if( mag(oldNormalContactForce_[cpI]) < smallFcn_ )
    //             {
    //                 activeFriction_[cpI] = false;
    //                 stick_[cpI] = false;
		    
    //                 oldRM_[cpI] =
    //                     beam.beamSegmentData
    //                     (
    //                         beam.solutionRMC(),
    //                         // beam.solutionRMC().oldTime(),
    //                         bI,
    //                         segI
    //                     )()[0];

    //                 torsionTheta_[cpI] = 0;
    //                 torsionDTheta_[cpI] = 0;

    //                 tangentialGap_[cpI] = 0;
    //                 tangentialGapOffset_[cpI] = 0;

    //                 axialTangentialGap_[cpI] = 0;
    //                 axialTangentialGapIncrement_[cpI] = 0;
    //                 axialTangentialGapOffset_[cpI] = 0;

    //                 transversalTangentialGap_[cpI] = 0;
    //                 transversalTangentialGapIncrement_[cpI] = 0;
    //                 transversalTangentialGapOffset_[cpI] = 0;

    //                 transversalFrictionalContactForce_ = vector::zero;
    //                 transversalFrictionalContactForceDerivativePerDTheta_ =
    //                     tensor::zero;
    //                 transversalFrictionalContactForceDerivativePerDW_ =
    //                     tensor::zero;

    //                 transversalFrictionalContactMoment_ = vector::zero;
    //                 transversalFrictionalContactMomentDerivativePerDTheta_ =
    //                     tensor::zero;
    //                 transversalFrictionalContactMomentDerivativePerDW_ =
    //                     tensor::zero;
    //             }
    //             else
    //             {
    //                 axialTangentialGap_[cpI] +=
    //                     axialTangentialGapIncrement_[cpI];

    //                 // axialTangentialGapIncrement_[cpI] = 0;

    //                 transversalTangentialGap_[cpI] +=
    //                     transversalTangentialGapIncrement_[cpI];
    //                 torsionTheta_[cpI] += torsionDTheta_[cpI];

    //                 // transversalTangentialGapIncrement_[cpI] = 0;

    //                 oldRM_[cpI] =
    //                     beam.beamSegmentData
    //                     (
    //                         beam.solutionRMC(),
    //                         // beam.solutionRMC().oldTime(),
    //                         bI,
    //                         segI
    //                     )()[0];

    //     	    // Axial tangential gap
    //                 axialTangentialGapOffset_[cpI] =
    //                     frictionModel.newGapOffset
    //                     (
    //                         axialTangentialGap_[cpI],
    //                         axialTangentialGapIncrement_[cpI],
    //                         mag(oldNormalContactForce_[cpI]),
    //                         axialTangentialGapOffset_[cpI]
    //                     );
    //                 axialTangentialGapIncrement_[cpI] = 0;

    //     	    // Transversal tangential gap
    //                 transversalTangentialGapOffset_[cpI] =
    //                     frictionModel.newGapOffset
    //                     (
    //                         transversalTangentialGap_[cpI],
    //                         transversalTangentialGapIncrement_[cpI],
    //                         mag(oldNormalContactForce_[cpI]),
    //                         transversalTangentialGapOffset_[cpI],
    //                         false
    //                     );
    //                 transversalTangentialGapIncrement_[cpI] = 0;
    //             }
    //         }
    //     }
    // }

    // const vectorField& cylCoord = pulleyContactCoordinates();

    // Tangential contact
    if (pulleyI != -1 && !frictionless)
    {
        tensor curRM =
            beam.beamSegmentData
            (
                beam.solutionRMC(),
                bI,
                segI
            )()[0];

        vector DW =
            beam.beamSegmentData
            (
                beam.solutionDW(),
                bI,
                segI
            )()[0];

        // vector W =
        //     beam.beamSegmentData
        //     (
        //         beam.solutionW(),
        //         bI,
        //         segI
        //     )()[0];

        // vector oldW =
        //     beam.beamSegmentData
        //     (
        //         beam.solutionW().oldTime(),
        //         bI,
        //         segI
        //     )()[0];

        // vector DW = W - oldW;

        // vector T =
        //     beam.beamSegmentData(beam.solutionCellTangent(), bI, segI)()[0];

        // vector T = splines[bI].dRdS()[segI];
        // T /= mag(T) + SMALL;

        vector T(1, 0, 0); // reference configuration
        T = (curRM & T); // current configuration
        T /= mag(T) + SMALL;

        // vector TT = (curRM.T() & T);
        // TT.y() = 0;
        // TT.z() = 0;
        // T = (curRM & TT);
        // T /= mag(T) + SMALL;

        // vector Tm = (curRM.T() & T);
        // Tm = (oldRM_[0] & Tm);
        // Tm += T;
        // Tm /= mag(Tm) + SMALL;

        // T0 = T;
        // vector T0(1, 0, 0);
        // T0 = (oldRM_[0] & T0);
        // T0 /= mag(T0) + SMALL;

        for (label cpI=0; cpI<2; cpI++)
        {
            if
            (
                activeFriction_[cpI]
             // && (mag(normalContactForce_[cpI]) > smallFcn_)
            )
            {
                if (true)
                // if (beam.iOuterCorr() > 0)
                {
                    vector tangentialGapIncrement = vector::zero;

                    if (mag(beam.U(bI)) > SMALL)
                    {
                        tangentialGapIncrement =
                           -beam.U(bI)*beam.runTime().deltaT().value()*T;

                        // Info << beam.U(bI) << endl;
                    }
                    else
                    {
                        tangentialGapIncrement = -DW;
                    }

                    // note - sign

                    // tangentialGapIncrement -=
                    //     normalDir_[cpI]*(normalDir_[cpI] & DW);

                    vector cylContactPointCoords =
                        pulleyContactCoordinates()[cpI];

                    // scalar rMin =
                    //     beam.conicalPulleys()[pulleyI]
                    //    .fullContactPointRadius(beam.R(bI));

                    // if (cylContactPointCoords.x()<rMin)
                    // {
                    //     cylContactPointCoords.x() = rMin;
                    // }

                    vector locPosition = 
                        beam.conicalPulleys()[pulleyI].position
                        (
                            cylContactPointCoords
                            // pulleyContactCoordinates()[cpI]
                        )
                      - beam.conicalPulleys()[pulleyI].origin();

                    scalar Omega = //0;
                        beam.conicalPulleys()[pulleyI]
                       .angularVelocity();
                       //  beam.conicalPulleys()[pulleyI]
                       // .angularVelocity(beam.R(bI));

                    // Info << "Omega: " << Omega << endl;

                    vector tangentialGapIncrementPulleyRot = // note + sign
                        (beam.conicalPulleys()[pulleyI].axis() ^ locPosition)
                        *Omega*beam.runTime().deltaT().value();
                    
                    // vector tangentialGapIncrement =
                    //     beam.conicalPulleys()[pulleyI].position
                    //     (
                    //         pulleyContactCoordinates()[cpI]
                    //     )
                    //   - beam.conicalPulleys()[pulleyI].position
                    //     (
                    //         oldPulleyContactCoordinates_[cpI]
                    //       // + vector(0, DPhi, 0)
                    //     );

                    if
                    (
                        // true
                           (mag(oldNormalContactForce_[0]) > smallFcn_)
                        && (mag(oldNormalContactForce_[1]) > smallFcn_)
                    )
                    {
                        axialTangentialGap_[cpI] = 0;
                        axialTangentialGapIncrement_[cpI] = 0;
                        axialFrictionalContactForce_[cpI] = vector::zero;
                        axialFrictionalContactForceDerivative_[cpI] =
                            tensor::zero;
                    }
                    else
                    {
                        scalar curAxialTangentialGap =
                            axialTangentialGap_[cpI]
                          + axialTangentialGapIncrement_[cpI];

                        scalar newAxialFrictionalContactForce =
                            frictionModel.contactForce
                            (
                                curAxialTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                axialTangentialGapOffset_[cpI]
                            );

                        axialFrictionalContactForce_[cpI] +=
                            frictionModel.axialRelaxationFactor()
                           *(
                               newAxialFrictionalContactForce
                             - (axialFrictionalContactForce_[cpI] & T)
                            )*T;

                        axialFrictionalContactForceDerivative_[cpI] =
                           -frictionModel.contactForceDerivative
                            (
                                curAxialTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                axialTangentialGapOffset_[cpI]
                            )*sqr(T);
                        
                        axialFrictionalContactForceDerivative_[cpI] +=
                            frictionModel.contactForceDerivative
                            (
                                curAxialTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                axialTangentialGapOffset_[cpI]
                            )
                           *Omega*beam.runTime().deltaT().value()
                           *(
                               (
                                   sqr(T)
                                   & spinTensor
                                   (
                                       beam.conicalPulleys()[pulleyI].axis()
                                       )
                                   )
                               & (I-sqr(normalDir_[cpI]))
                               );
                        
                        axialFrictionalContactForceDerivative_[cpI] +=
                            frictionModel.contactForceDerivativeOverFcn
                            (
                                curAxialTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                axialTangentialGapOffset_[cpI],
                                false
                            )
                           *normalModel.contactForceDerivative
                            (
                                normalGap_[cpI],
                                0 //normalContactForceOffset()
                            )
                           *(T*normalDir_[cpI]);
                    }

                    // Transversal fricion force

                    vector R0 = -normalDir_[cpI];
                    R0 -= T*(T&R0);
                    R0 /= mag(R0) + SMALL;
                    // transversal friction force postive direction
                    vector u0 = (T ^ R0);

                    scalar curTransversalTangentialGap =
                        transversalTangentialGap_[cpI]
                      + transversalTangentialGapIncrement_[cpI];

                    // force
                    if
                    (
                        // true
                        (
                            (mag(oldNormalContactForce_[0]) > smallFcn_)
                         && (mag(oldNormalContactForce_[1]) > smallFcn_)
                        )
                     && frictionlessFullContact
                    )
                    {
                        transversalTangentialGap_[cpI] = 0;
                        transversalTangentialGapIncrement_[cpI] = 0;
                        transversalTangentialGapOffset_[cpI] = 0;

                        transversalFrictionalContactForceDerivativePerDW_[cpI] =
                            tensor::zero;
                        transversalFrictionalContactForceDerivativePerDTheta_[cpI] =
                            tensor::zero;
                        transversalFrictionalContactForce_[cpI] =
                            vector::zero;
                    }
                    else
                    {
                        // scalar newTransversalFrictionalContactForce =
                        //     frictionModel.contactForce
                        //     (
                        //         curTransversalTangentialGap,
                        //         mag(normalContactForce_[cpI]) + smallFcn_,
                        //         transversalTangentialGapOffset_[cpI],
                        //         false
                        //     );
                    
                        transversalFrictionalContactForce_[cpI] =
                            frictionModel.contactForce
                            (
                                curTransversalTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                transversalTangentialGapOffset_[cpI],
                                false
                            )*u0;

                        // transversalFrictionalContactForce_[cpI] +=
                        //     frictionModel.transversalRelaxationFactor()
                        //    *(
                        //         newTransversalFrictionalContactForce
                        //       - (transversalFrictionalContactForce_[cpI] & u0)
                        //     )*u0;
                        
                        transversalFrictionalContactForceDerivativePerDTheta_[cpI] = 
                           -frictionModel.contactForceDerivative
                            (
                                curTransversalTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                transversalTangentialGapOffset_[cpI],
                                false
                                )*beam.R(bI)*(u0*T);
                        
                        transversalFrictionalContactForceDerivativePerDW_[cpI] =
                            frictionModel.contactForceDerivativeOverFcn
                            (
                                curTransversalTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                transversalTangentialGapOffset_[cpI],
                                false
                                )
                            *normalModel.contactForceDerivative
                            (
                                normalGap_[cpI],
                                0 //normalContactForceOffset()
                                )
                            *(u0*normalDir_[cpI]);
                        
                        // transversalFrictionalContactForceDerivativePerDW_[cpI] +=
                        //     frictionModel.contactForceDerivative
                        //     (
                        //         curTransversalTangentialGap,
                        //         mag(oldNormalContactForce_[cpI]) + smallFcn_,
                        //         transversalTangentialGapOffset_[cpI],
                        //         false
                        //     )
                        //    *Omega*beam.runTime().deltaT().value()
                        //    *(
                        //         (
                        //             sqr(u0)
                        //           & spinTensor
                        //             (
                        //                 beam.conicalPulleys()[pulleyI].axis()
                        //             )
                        //         )
                        //       & (I-sqr(normalDir_[cpI]))
                        //     );

                    
                        // transversalFrictionalContactForceDerivativePerDW_[cpI] =
                        //     tensor::zero;
                        // transversalFrictionalContactForceDerivativePerDTheta_[cpI] =
                        //     tensor::zero;
                    }
                
                    // moment
                    if
                    (
                        // true
                        (
                            (mag(oldNormalContactForce_[0]) > smallFcn_)
                         && (mag(oldNormalContactForce_[1]) > smallFcn_)
                        )
                     && frictionlessFullContact
                    )
                    {
                        transversalTangentialGap_[cpI] = 0;
                        transversalTangentialGapIncrement_[cpI] = 0;
                        transversalTangentialGapOffset_[cpI] = 0;

                        transversalFrictionalContactMomentDerivativePerDW_[cpI] =
                            tensor::zero;
                        transversalFrictionalContactMomentDerivativePerDTheta_[cpI] =
                            tensor::zero;
                        transversalFrictionalContactMoment_[cpI] =
                            vector::zero;
                    }
                    else
                    {
                        // scalar newTransversalFrictionalContactMoment = 
                        //     frictionModel.contactForce
                        //     (
                        //         curTransversalTangentialGap,
                        //         mag(normalContactForce_[cpI]) + smallFcn_,
                        //         transversalTangentialGapOffset_[cpI],
                        //         false
                        //     )*beam.R(bI);

                        transversalFrictionalContactMoment_[cpI] = 
                            frictionModel.contactForce
                            (
                                curTransversalTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                transversalTangentialGapOffset_[cpI],
                                false
                            )*beam.R(bI)*T;

                        // transversalFrictionalContactMoment_[cpI] +=
                        //     frictionModel.transversalRelaxationFactor()
                        //    *(
                        //         newTransversalFrictionalContactMoment
                        //       - (transversalFrictionalContactMoment_[cpI] & T)
                        //     )*T;
                        
                        transversalFrictionalContactMomentDerivativePerDTheta_[cpI] =
                           -frictionModel.contactForceDerivative
                            (
                                curTransversalTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                transversalTangentialGapOffset_[cpI],
                                false
                            )*sqr(beam.R(bI))*(T*T);

                        transversalFrictionalContactMomentDerivativePerDW_[cpI] =
                            frictionModel.contactForceDerivativeOverFcn
                            (
                                curTransversalTangentialGap,
                                mag(normalContactForce_[cpI]) + smallFcn_,
                                transversalTangentialGapOffset_[cpI],
                                false
                            )
                           *normalModel.contactForceDerivative
                            (
                                normalGap_[cpI],
                                0 //normalContactForceOffset()
                            )
                           *beam.R(bI)*(T*normalDir_[cpI]);

                        // transversalFrictionalContactMomentDerivativePerDW_[cpI] +=
                        //     frictionModel.contactForceDerivative
                        //     (
                        //         curTransversalTangentialGap,
                        //         mag(normalContactForce_[cpI]) + smallFcn_,
                        //         transversalTangentialGapOffset_[cpI],
                        //         false
                        //     )
                        //    *Omega*beam.runTime().deltaT().value()
                        //    *beam.R(bI)
                        //    *(
                        //         (
                        //             (T*u0)
                        //           & spinTensor
                        //             (
                        //                 beam.conicalPulleys()[pulleyI].axis()
                        //             )
                        //         )
                        //       & (I-sqr(normalDir_[cpI]))
                        //     );
                    
                        // transversalFrictionalContactMomentDerivativePerDW_[cpI] =
                        //     tensor::zero;
                        // transversalFrictionalContactMomentDerivativePerDTheta_[cpI] =
                        //     tensor::zero;
                    }     
                }
                else
                {
                    // Info << "Predictor: " << beam.iOuterCorr() << endl;
                }
            }
            else
            {
                // frictionalContactMoment_[cpI] = vector::zero;
                // frictionalContactMomentDerivative_[cpI] = tensor::zero;
                // frictionalContactMomentDerivativeOverFcn_[cpI] = tensor::zero;
                // frictionalContactForce_[cpI] = vector::zero;
                // frictionalContactForceDerivative_[cpI] = tensor::zero;
                // frictionalContactForceDerivativeOverFcn_[cpI] = tensor::zero;

                axialFrictionalContactForce_[cpI] = vector::zero;
                axialFrictionalContactForceDerivative_[cpI] = tensor::zero;

                transversalFrictionalContactForce_[cpI] = vector::zero;
                transversalFrictionalContactForceDerivativePerDTheta_[cpI] =
                    tensor::zero;
                transversalFrictionalContactForceDerivativePerDW_[cpI] =
                    tensor::zero;
                
                transversalFrictionalContactMoment_[cpI] = vector::zero;
                transversalFrictionalContactMomentDerivativePerDTheta_[cpI] =
                    tensor::zero;
                transversalFrictionalContactMomentDerivativePerDW_[cpI] =
                    tensor::zero;
            }
        }
    }
}
// ************************************************************************* //
