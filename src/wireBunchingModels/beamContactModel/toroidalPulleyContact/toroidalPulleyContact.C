/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield160
         | foam-extend: Open Source CFD
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

#include "toroidalPulleyContact.H"
#include "beamModel.H"
// #include "scalarMatrices.H"
#include "pseudoVector.H"
#include "spinTensor.H"
#include "lagrangeMultipliers.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::toroidalPulleyContact::smallFcn_ = 1e-3;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::toroidalPulleyContact::toroidalPulleyContact()
:
    firstBeam_(-1),
    firstBeamSegment_(-1),
    firstBeamZeta_(0),
    pulleyIndex_(-1),
    pulleyContactCoordinate_(vector::zero),
    delta_(GREAT),
    normalDir_(vector::zero),
    contactAngle_(0),
    normalContactForce_(vector::zero),
    oldNormalContactForce_(vector::zero),
    normalContactForceDerivative_(tensor::zero),
    normalContactForceDirectionDerivative_(tensor::zero),
    normalGap_(GREAT),
    oldNormalGap_(GREAT),
    frictionalContactMoment_(vector::zero),
    frictionalContactForce_(vector::zero),
    axialFrictionalContactForce_(vector::zero),
    frictionalContactMomentDerivative_(tensor::zero),
    frictionalContactMomentDerivativeOverFcn_(tensor::zero),
    frictionalContactForceDerivative_(tensor::zero),
    frictionalContactForceDerivativeOverFcn_(tensor::zero),
    oldTorsionTheta_(0),
    torsionDTheta_(0),
    oldRM_(tensor::I),
    tangentialGap_(0),
    tangentialGapOffset_(0),
    activeFriction_(false),
    stick_(true),
    unloading_(false),
    curTimeIndex_(-1)
{}

Foam::toroidalPulleyContact::toroidalPulleyContact
(
    const label fb,
    const label fbSeg,
    const scalar fbZeta,
    const label pulleyIndex,
    const vector pcc,
    const label timeIndex
)
:
    firstBeam_(fb),
    firstBeamSegment_(fbSeg),
    firstBeamZeta_(fbZeta),
    pulleyIndex_(pulleyIndex),
    pulleyContactCoordinate_(pcc),
    delta_(GREAT),
    normalDir_(vector::zero),
    contactAngle_(0),
    normalContactForce_(vector::zero),
    oldNormalContactForce_(vector::zero),
    normalContactForceDerivative_(tensor::zero),
    normalContactForceDirectionDerivative_(tensor::zero),
    normalGap_(GREAT),
    oldNormalGap_(GREAT),
    frictionalContactMoment_(vector::zero),
    frictionalContactForce_(vector::zero),
    axialFrictionalContactForce_(vector::zero),
    frictionalContactMomentDerivative_(tensor::zero),
    frictionalContactMomentDerivativeOverFcn_(tensor::zero),
    frictionalContactForceDerivative_(tensor::zero),
    frictionalContactForceDerivativeOverFcn_(tensor::zero),
    oldTorsionTheta_(0),
    torsionDTheta_(0),
    // oldTorsionDTheta_(2, 0),
    oldRM_(tensor::I),
    tangentialGap_(0),
    tangentialGapOffset_(0),
    activeFriction_(false),
    stick_(true),
    unloading_(false),
    curTimeIndex_(timeIndex)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::toroidalPulleyContact::set
(
    const label fb,
    const label fbSeg,
    const scalar fbZeta,
    const label pulleyIndex,
    const vector pulleyContactCoordinate,
    const label timeIndex
)
{
    firstBeam_ = fb;
    firstBeamSegment_ = fbSeg;
    firstBeamZeta_ = fbZeta;

    pulleyIndex_ = pulleyIndex;
    pulleyContactCoordinate_ = pulleyContactCoordinate;
}


void Foam::toroidalPulleyContact::updateForce
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

    if (curTimeIndex_ < beam.runTime().timeIndex())
    {
        curTimeIndex_ = beam.runTime().timeIndex();

        oldNormalContactForce_ = normalContactForce_;
        oldNormalGap_ = normalGap_;

        // // Activate/deactivate friction
	// for (label cpI=0; cpI<2; cpI++)
	// {
        //     scalar maxFrictionForce =
        //         frictionModel.maxFrictionForce
        //         (
        //             mag(oldNormalContactForce_[cpI])
        //         );

        //     if (!activeFriction_[cpI])
        //     {
        //         if
        //         (
        //             mag(oldNormalContactForce_[cpI]) > smallFcn_
        //          && maxFrictionForce > SMALL
        //         )
        //         {
        //             activeFriction_[cpI] = true;
        //             stick_[cpI] = true;

        //             // Contact inception: set initial RM
        //             oldRM_ =
        //                 beam.beamSegmentData
        //                 (
        //                     beam.solutionRMC(),
        //                     bI,
        //                     segI
        //                 )()[0];

        //             oldTorsionTheta_[cpI] = 0;
        //             torsionDTheta_[cpI] = 0;
        //             tangentialGapOffset_[cpI] = 0;
        //         }
        //     }
        //     else
        //     {
        //         if
        //         (
        //             mag(oldNormalContactForce_[cpI]) < smallFcn_
        //         )
        //         {
        //             activeFriction_[cpI] = false;
        //             stick_[cpI] = true;
        //             oldRM_[cpI] = tensor::zero;
        //             oldTorsionTheta_[cpI] = 0;
        //             torsionDTheta_[cpI] = 0;
        //             tangentialGapOffset_[cpI] = 0;
        //         }
        //         else
        //         {
        //             oldRM_ =
        //                 beam.beamSegmentData
        //                 (
        //                     beam.solutionRMC(),
        //                     bI,
        //                     segI
        //                 )()[0];

        //             oldTorsionTheta_[cpI] += torsionDTheta_[cpI];

        //             scalar oldGap = beam.R(bI)*oldTorsionTheta_[cpI];

        //             const bool stick =
        //                 frictionModel.stick
        //                 (
        //                     oldGap,
        //                     mag(oldNormalContactForce_[cpI]),
        //                     tangentialGapOffset_[cpI]
        //                 );

        //             if (!stick) // slip
        //             {
        //                 scalar DGap = beam.R(bI)*torsionDTheta_[cpI];
        //                 DGap *= oldGap/(mag(oldGap) + SMALL);

        //                 if (DGap < -SMALL)
        //                 {
        //                     unloading_[cpI] = true;
        //                 }
        //                 else
        //                 {
        //                     unloading_[cpI] = false;
        //                 }
        //             }

        //             tangentialGapOffset_[cpI] =
        //                 frictionModel.newGapOffset
        //                 (
        //                     oldGap,
        //                     mag(oldNormalContactForce_[cpI]),
        //                     tangentialGapOffset_[cpI]
        //                 ); 
        //         }
        //     }
        // }
    }

    label pulleyI = pulleyIndex();
    const vector& cylCoord = pulleyContactCoordinate();

    vector cylOrigin = cylCoord;
    cylOrigin.x() = 0;
    
    scalar dist = GREAT;
    vector n = vector::one;

    if (pulleyI != -1)
    {
        // vector T = splines[bI].dRdS()[segI];
        // T /= mag(T) + SMALL;
        
	{
            vector curContactPoint = 
                splines[bI].position(segI, zeta);

	    vector curNeiContactPoint = 
	        beam.toroidalPulleys()[pulleyI].position(cylCoord);

	    vector origin = 
	        beam.toroidalPulleys()[pulleyI].position(cylOrigin);

            n = curNeiContactPoint - origin;
	    n /= mag(n) + SMALL;

	    dist = ((curContactPoint - curNeiContactPoint) & n);

	    // n = curContactPoint - curNeiContactPoint;
	    // n /= mag(n) + SMALL;

	    scalar gap = dist - beam.R(bI);

	    normalGap_ = gap;

	    delta() = cylCoord.x();

            // // Subtract axial component
            // n = n - T*(T & n);
            // n /= mag(n) + SMALL;

            normalDir_ = n;
            
            if (!isA<lagrangeMultipliers>(normalModel))
            {
                normalContactForce() =
                    normalModel.contactForce
                    (
                        gap,
                        0 //normalContactForceOffset()
                    )*n;


                normalContactForceDerivative() =
                    normalModel.contactForceDerivative
                    (
                        gap,
                        0 //normalContactForceOffset()
                    )*sqr(n);

                // if (pulleyI == 5)
                // {
                //     if ( mag(normalContactForce()) > SMALL)
                //     {
                //         scalar fcder =
                //             normalModel.contactForceDerivative
                //             (
                //                 gap,
                //                 0 //normalContactForceOffset()
                //             );
                        
                //         Info << segI << " " << gap << " "
                //              << normalContactForce() << " "
                //              << fcder  << endl;
                //         Info << cylCoord << endl;
                //     }
                // }

                // vector DrDz =
                //     beam.toroidalPulleys()[pulleyI].DrDz(cylCoord)
                //    *beam.toroidalPulleys()[pulleyI].axis();

                // normalContactForceDerivative() -=
                //     normalModel.contactForceDerivative
                //     (
                //         gap,
                //         0
                //     )*(n*DrDz);

                vector t = vector(0, 0, 1);

                normalContactForceDirectionDerivative() = //tensor::zero;
                    mag(normalContactForce())
                   *((I-sqr(t))&(I-sqr(n)))/delta();
                
                if (agumentedLagrangian)
                {
                    // if ((gap < -1e-4*beam.R().value()) || (gap > 0))
                    // {
                    // normalContactForceOffset() +=
		    // (    
                    //     mag(normalContactForce())
		    //   - normalContactForceOffset()
		    // );
                    // }
                }
            }
	}
    }
    else
    {
        normalContactForce() = vector::zero;
        normalContactForceDerivative() = tensor::zero;
        normalContactForceDirectionDerivative() = tensor::zero;
    }
}

// ************************************************************************* //
