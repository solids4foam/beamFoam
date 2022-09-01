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

#include "lineContact.H"
#include "beamModel.H"
#include "scalarMatrices.H"
#include "pseudoVector.H"
#include "spinTensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineContact::lineContact()
:
    firstBeam_(-1),
    firstBeamSegment_(-1),
    firstBeamOldSegment_(-1),
    firstBeamZeta_(0),
    firstBeamOldZeta_(0),
    contactDirection_(vector::zero),
    oldContactDirection_(vector::zero),
    secondBeam_(-1),
    secondBeamSegment_(-1),
    secondBeamLowerSegment_(-1),
    secondBeamUpperSegment_(-1),
    secondBeamWeight_(0),
    secondBeamOldSegment_(-1),
    secondBeamZeta_(0),
    secondBeamOldZeta_(0),
    delta_(GREAT),
    contactAngle_(90),
    normalContactForce_(vector::zero),
    prevNormalContactForce_(vector::zero),
    normalContactForceOffset_(0),
    normalContactForceDerivative_(tensor::zero),
    axialContactForce_(vector::zero),
    axialContactForceDerivative_(tensor::zero),
    circumferentialContactForce_(vector::zero),
    circumferentialContactMoment_(vector::zero),
    circumferentialContactForceDerivative_(tensor::zero),
    circumferentialContactForceThetaDerivative_(tensor::zero),
    circumferentialContactMomentDerivative_(tensor::zero),
    normalGap_(GREAT),
    normalGapOffset_(0),
    axialGap_(0),
    oldAxialGap_(0),
    circumferentialGap_(0),
    oldCircumferentialGap_(0),
    axialSlipGap_(0),
    oldAxialSlipGap_(0),
    circumferentialSlipGap_(0),
    oldCircumferentialSlipGap_(0),
    curTimeIndex_(-1)
{}

Foam::lineContact::lineContact
(
    const label fb,
    const label fbSeg,
    const scalar fbZeta,
    const label sb,
    const label sbSeg,
    const scalar sbZeta,
    const label timeIndex
)
:
    firstBeam_(fb),
    firstBeamSegment_(fbSeg),
    firstBeamOldSegment_(-1),
    firstBeamZeta_(fbZeta),
    firstBeamOldZeta_(0),
    contactDirection_(vector::zero),
    oldContactDirection_(vector::zero),
    secondBeam_(sb),
    secondBeamSegment_(sbSeg),
    secondBeamLowerSegment_(-1),
    secondBeamUpperSegment_(-1),
    secondBeamWeight_(0),
    secondBeamOldSegment_(-1),
    secondBeamZeta_(sbZeta),
    secondBeamOldZeta_(0),
    delta_(GREAT),
    contactAngle_(90),
    normalContactForce_(vector::zero),
    prevNormalContactForce_(vector::zero),
    normalContactForceOffset_(0),
    normalContactForceDerivative_(tensor::zero),
    axialContactForce_(vector::zero),
    axialContactForceDerivative_(tensor::zero),
    circumferentialContactForce_(vector::zero),
    circumferentialContactMoment_(vector::zero),
    circumferentialContactForceDerivative_(tensor::zero),    
    circumferentialContactForceThetaDerivative_(tensor::zero),
    circumferentialContactMomentDerivative_(tensor::zero),
    normalGap_(GREAT),
    normalGapOffset_(0),
    axialGap_(0),
    oldAxialGap_(0),
    circumferentialGap_(0),
    oldCircumferentialGap_(0),
    axialSlipGap_(0),
    oldAxialSlipGap_(0),
    circumferentialSlipGap_(0),
    oldCircumferentialSlipGap_(0),
    curTimeIndex_(timeIndex)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lineContact::set
(
    const label fb,
    const label fbSeg,
    const scalar fbZeta,
    const label sb,
    const label sbSeg,
    const scalar sbZeta,
    const label timeIndex
)
{
	// Info << "Inside set() of lineContact.C" << endl;
    if (curTimeIndex_ < timeIndex)
    {
        if (secondBeamOldSegment_ != -1) //There is contact
		{
			// Info << "There is contact" << endl;
			
            firstBeamOldSegment_ = firstBeamSegment_;
			firstBeamOldZeta_ = firstBeamZeta_;

			oldContactDirection_ = contactDirection_;

			secondBeamOldSegment_ = secondBeamSegment_;
			secondBeamOldZeta_ = secondBeamZeta_;

			oldAxialGap_ = axialGap_;
			oldCircumferentialGap_ = circumferentialGap_;

			oldAxialSlipGap_ = axialSlipGap_;
			oldCircumferentialSlipGap_ = circumferentialSlipGap_;
		}
		
		// Info << "No contact" << endl;
	
        curTimeIndex_  = timeIndex;
    }

    firstBeam_ = fb;
    firstBeamSegment_ = fbSeg;
    firstBeamZeta_ = fbZeta;

    secondBeam_ = sb;
    secondBeamSegment_ = sbSeg;
    secondBeamZeta_ = sbZeta;
}


void Foam::lineContact::updateForce
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines,
    const normalContactModel& normalModel,
    const frictionContactModel& frictionModel,
    const scalar lowerContactAngleLimit,
    const scalar upperContactAngleLimit,
    const bool agumentedLagrangian
)
{
    // bool debug = false;
    // label debugSeg = 99;

    label bI = firstBeam();
    label segI = firstBeamSegment();
    scalar zeta = firstBeamZeta();

    label nbI = secondBeam();
    label neiSegI = secondBeamSegment();
    scalar neiZeta = secondBeamZeta();

    scalar dist = GREAT;
    vector n = vector::one;

    
    // Determine upper and lower segments (cells)
    if (neiSegI != -1)
    {
        const label nSegments = beam.contact().splines()[nbI].nSegments();

        if
        (
            (neiSegI == 0 && neiZeta <= 0)
         || (neiSegI == (nSegments-1) && neiZeta >= 0)
        )
        {
            secondBeamLowerSegment_ = neiSegI;
            secondBeamUpperSegment_ = neiSegI;
        }
        else if (neiZeta > 0)
        {
            secondBeamLowerSegment_ = neiSegI;
            secondBeamUpperSegment_ = neiSegI+1;
        }
        else
        {
            secondBeamLowerSegment_ = neiSegI-1;
            secondBeamUpperSegment_ = neiSegI;
        }

        if (secondBeamLowerSegment_ == secondBeamUpperSegment_)
        {
            secondBeamWeight_ = 1.0;
        }
        else
        {
            scalar l =
                0.5*beam.contact().splines()[nbI]
               .segLength(secondBeamLowerSegment_)
              + 0.5*beam.contact().splines()[nbI]
               .segLength(secondBeamUpperSegment_);

            scalar lx =
                mag(neiZeta)
               *beam.contact().splines()[nbI]
               .segLength(neiSegI)/2;

            if (neiSegI == secondBeamLowerSegment_)
            {
                secondBeamWeight_ = (l-lx)/l;
            }
            else
            {
                secondBeamWeight_ = 1.0 - (l-lx)/l;
            }            
        }
    }


    if (neiSegI != -1)
    {
        // Calculating contact angle
        vector tan0 = splines[bI].firstDerivative(segI, 0);
        
	tan0 /= mag(tan0) + SMALL;

	vector tan1 =
		splines[nbI].firstDerivative(neiSegI, neiZeta);
			
	tan1 /= mag(tan1) + SMALL;

	contactAngle_ = ::acos(mag(tan0 & tan1))*180.0/M_PI;
	
	//   Info << "contact angle lineContact.C file " 
	     //   << ::acos(mag(tan0 & tan1))*180.0/M_PI 
	     //   << endl;
		
	// if (contactAngle_ < upperContactAngleLimit)
	// {
	    // SB: why is the contact point on main beam
	    // not calculated at the centre of CV
            vector curContactPoint = 
                splines[bI].position(segI, zeta);
		
	    // Info << "zeta " << zeta << endl;

	    vector curNeiContactPoint = 
		splines[nbI].position(neiSegI, neiZeta);

	    dist = mag(curContactPoint - curNeiContactPoint);

	    n = curContactPoint - curNeiContactPoint;
	    n /= mag(n) + SMALL;

            // Calculate contact gap (penetration)
	    scalar gap = dist - beam.R(bI) - beam.R(nbI);
	    

            // Determine conection line ref angle
            if (!Pstream::parRun())
            {
                tensor curRM =
                    beam.beamSegmentData
                    (
                        beam.solutionRMC(),
                        bI,
                        segI
                    )()[0];
                
                vector dHat =
                    curNeiContactPoint - curContactPoint;
                dHat /= mag(dHat) + SMALL;

                vector dRefHat = (curRM.T() & dHat);
                dRefHat.x() = 0;
                dRefHat /= mag(dRefHat);

                scalar theta =
                    beam.crossSections()[bI]
                   .polarAngle
                    (
                        dRefHat.z(),
                        dRefHat.y()
                    );

                // Determine conection line ref angel
                tensor curNeiRM =
                    beam.beamSegmentData
                    (
                        beam.solutionRMC(),
                        nbI,
                        neiSegI
                    )()[0];

                vector neiDHat =
                    curContactPoint - curNeiContactPoint;
                neiDHat /= mag(neiDHat) + SMALL;

                vector neiDRefHat = (curNeiRM.T() & neiDHat);
                neiDRefHat.x() = 0;
                neiDRefHat /= mag(neiDRefHat);

                scalar newTheta =
                    beam.crossSections()[nbI]
                   .polarAngle
                    (
                        neiDRefHat.z(),
                        neiDRefHat.y()
                    );

                // Calculate contact gap (penetration)
                scalar gap =
                    dist
                  - beam.crossSections()[bI].radius(theta)
                  - beam.crossSections()[nbI].radius(newTheta);
		  
	      
            }

            // if (gap < 0)
            // {
            //     Info << theta << ", " << newTheta << endl;
            // }

	    normalGap_ = gap;

	    delta() = dist;

	    prevNormalContactForce_ = normalContactForce();
	    
	    scalar lineWeight = 
		allAngleLineContactWeight
		(
		    beam,
		    splines,
		    lowerContactAngleLimit,
		    upperContactAngleLimit
		);
		
	    // Info << "lineWeight " << lineWeight << endl;

	    normalContactForce() = lineWeight*
		normalModel.contactForce
		    (
			gap,
			normalContactForceOffset(),
			normalGapOffset()
		    )*n;
		    
	    //Info << "contact force " << mag(normalContactForce()) << endl;

	    if (agumentedLagrangian)
	    {
		if (mag(prevNormalContactForce_) > SMALL)
	        // if ((gap < -1e-4*(beam.R(bI)+beam.R(nbI))/2) || (gap > 0))
		{
		    normalContactForceOffset() +=
			0.5*
			(    
				mag(normalContactForce())
				- normalContactForceOffset()
			);

		    normalContactForce() =
			normalModel.contactForce
			(
				gap,
				normalContactForceOffset()
			)*n;
		}
	    }

	    normalContactForceDerivative() = lineWeight*
		normalModel.contactForceDerivative
		    (
			gap,
			normalContactForceOffset(),
			normalGapOffset()
		    )*sqr(n);
	// }
	//else
	//{
   	    // // Info << "contactAngle: " << contactAngle_ << endl;

		// normalContactForce() = vector::zero;
		// normalContactForceDerivative() = tensor::zero;
	//}
    }
    else
    {
        normalContactForce() = vector::zero;
        normalContactForceDerivative() = tensor::zero;
    }

    // Tangential force
    if
    (
        (mag(normalContactForce_) > SMALL)
		&& (secondBeamOldSegment_ == -1)
		&& false
    )
    {
        firstBeamOldSegment_ = firstBeamSegment();
        firstBeamOldZeta_ = firstBeamZeta();

        secondBeamOldSegment_ = secondBeamSegment();
        secondBeamOldZeta_ = secondBeamZeta();

		oldContactDirection_ = contactDirection_;
    }

    // Here this entire loop is not executed. 
    // Hence, the axial and circumferential contact forces are zero
    // See the else loop at around line no 872
	if
	(
		(mag(normalContactForce_) > SMALL)&& false
	)
	{	
		scalar maxFrictionForce =
            frictionModel.maxFrictionForce
            (
                mag(normalContactForce())
            );

		// Force in axial direction
        label oldNeiSegI = secondBeamOldSegment_;
        scalar oldNeiZeta = secondBeamOldZeta_;

		vector oldNeiPoint =
			splines[nbI].position(oldNeiSegI, oldNeiZeta);
		labelScalar np =
			splines[bI].findNearestPoint
			(
				segI,
				oldNeiPoint,
				2
			);

		label oldSegI = np.first();
		scalar oldZeta = np.second();

        axialGap_ =
            splines[bI].length
            (
                oldSegI,
                oldZeta
            )
          - splines[bI].length
            (
                segI,
                zeta // zero
            );
		axialGap_ += oldAxialGap_;

        vector T = splines[bI].paramFirstDerivative(segI, zeta);
        scalar J = mag(T);
        T /= J;

        vector neiT = splines[nbI].paramFirstDerivative(neiSegI, neiZeta);
        scalar neiJ = mag(neiT);
        neiT /= neiJ;

        scalar trialElasticAxialGap = axialGap_ - oldAxialSlipGap_;

        vector trialAxialContactForce = 
            frictionModel.contactForce
            (
                mag(trialElasticAxialGap)
            )*T*sign(trialElasticAxialGap);

        if (mag(trialAxialContactForce) > maxFrictionForce)
        {
            axialContactForce_ =
                maxFrictionForce*T*sign(trialElasticAxialGap);

			scalar forceCorrection =
			(
				(sign(trialElasticAxialGap)*T)
				& (trialAxialContactForce - axialContactForce_)
			);

            axialSlipGap_ =
                oldAxialSlipGap_
              + frictionModel.slipGapCorrection
                (
                    mag(trialElasticAxialGap),
                    forceCorrection
                )*sign(trialElasticAxialGap);

            trialElasticAxialGap = axialGap_ - axialSlipGap_;

            trialElasticAxialGap = 0; // Zero contact force derivative
        }
        else if (mag(trialAxialContactForce) < SMALL)
        {
            axialContactForce_ =
                maxFrictionForce*T*sign(trialElasticAxialGap);

            axialSlipGap_ = axialGap_;
            trialElasticAxialGap = axialGap_ - axialSlipGap_;
        }
        else
        {
            axialContactForce_ = trialAxialContactForce;
            axialSlipGap_ = oldAxialSlipGap_;
        }

	    // if
	    // (
	    //     (bI == 0)
	    //  && (segI == debugSeg)
	    //  && debug
            // )
	    // {
	    //     Info << "a1" << endl;
	    //     Info << "axialGap: " << axialGap_ << endl;
	    //     Info << "axialSlipGap: " << axialSlipGap_ << endl;
	    //     Info << "oldAxialSlipGap: " << oldAxialSlipGap_ << endl;
	    //     Info << "normalContactForce: "
	    //          << mag(normalContactForce_) << endl;
	    // }
	    
		// axialContactForce_ = vector::zero;

	axialContactForceDerivative_ = //tensor::zero;
		sign(trialElasticAxialGap)*(-1)
			*frictionModel.contactForceDerivative
            (
                mag(trialElasticAxialGap)
            )*(T*T);

		// Force in circumferential direction

	// Circumferetial gap due to relative twisting of beams
	contactDirection_ = -normalContactForce_;
	contactDirection_ /= mag(contactDirection_);

	vector d = contactDirection_;
	d -= T*(T&d);
	d /= mag(d) + SMALL;

	vector d0 = oldContactDirection_;
	d0 -= T*(T&d0);
	d0 /= mag(d0) + SMALL;

	scalar twistDAlpha = ::asin(mag(d0 ^ d));
	scalar twistSgn = -sign(T & (d0 ^ d));
	scalar twistCircumGap =
	    twistSgn*(2*twistDAlpha*beam.R(bI)); // ZT: check this (bI or nbI)


	// Circumferential gap due to relative torsion of beams
	const surfaceTensorField& Lambda =
	    beam.solutionLambda();

	tensorField segLambda =
	    beam.beamPointData(Lambda, bI, segI);
	tensor dLambda = (segLambda[0].T() & segLambda[1]);
	vector dTheta = pseudoVector(dLambda);

	tensor LambdaC =
	    (segLambda[0] & rotationMatrix(0.5*dTheta));

	vector R = d*beam.R(bI); // ZT: check this (bI or nbI)
	vector R0 = (LambdaC.T() & R);

	R -= T*(T&R);
	R /= mag(R) + SMALL;

	R0 -= T*(T&R0);
	R0 /= mag(R0) + SMALL;

	scalar DAlpha = ::asin(mag(R0 ^ R));
	scalar sgn = sign(T & (R0 ^ R));

	tensorField neiSegLambda =
	    beam.beamPointData(Lambda, nbI, neiSegI);
	tensor neiDLambda = (neiSegLambda[0].T() & neiSegLambda[1]);
	vector neiDTheta = pseudoVector(neiDLambda);

	tensor neiLambdaC =
	    (neiSegLambda[0] & rotationMatrix((neiZeta+1)*neiDTheta/2));

	vector neiR = -d*beam.R(nbI); // ZT: check this (bI or nbI)
	vector neiR0 = (neiLambdaC.T() & neiR);

	neiR -= neiT*(neiT&neiR);
	neiR /= mag(neiR) + SMALL;

	neiR0 -= neiT*(neiT&neiR0);
	neiR0 /= mag(neiR0) + SMALL;

	scalar neiDAlpha = ::asin(mag(neiR0 ^ neiR));
	scalar neiSgn = sign(neiT & (neiR0 ^ neiR));

	scalar torsionCircumGap =
	    sgn*(DAlpha*beam.R(nbI)) // ZT: check this (bI or nbI)
  	  + neiSgn*(neiDAlpha*beam.R(nbI))*::cos(contactAngle_*M_PI/180);

	// torsionCircumGap = 0;
	// twistCircumGap = 0;

	circumferentialGap_ =
	    torsionCircumGap
	  + twistCircumGap
	  + oldCircumferentialGap_;

        scalar trialElasticCircumGap =
	    circumferentialGap_
	  - oldCircumferentialSlipGap_;

	vector circumDir = (T ^ n);
	// vector circumDir = (T ^ d);

        vector trialCircumContactForce =
            frictionModel.contactForce
            (
                mag(trialElasticCircumGap)
	    )*circumDir*sign(trialElasticCircumGap);
 
	vector prevCircumferentialContactForce =
	    circumferentialContactForce_;

        if (mag(trialCircumContactForce) > maxFrictionForce)
        {
            circumferentialContactForce_ =
                maxFrictionForce*circumDir
	       *sign(trialElasticCircumGap);

	    scalar forceCorrection =
	    (
	        (circumDir*sign(trialElasticCircumGap))
	      & (trialCircumContactForce - circumferentialContactForce_)
	    );

            circumferentialSlipGap_ =
                oldCircumferentialSlipGap_
              + frictionModel.slipGapCorrection
                (
                    mag(trialElasticCircumGap),
		    forceCorrection
                )*sign(trialElasticCircumGap);

            trialElasticCircumGap =
	        circumferentialGap_
	      - circumferentialSlipGap_;

	    // if
	    // (
	    //     (bI == 0)
	    //  && (segI == debugSeg)
	    //  && debug
            // )
	    // {
	    //     Info << "x1" << endl;
	    //     Info << circumferentialGap_ << endl;
	    //     Info << circumferentialSlipGap_ << endl;
	    //     Info << oldCircumferentialSlipGap_ << endl;
	    //     Info << "trialElasticCircumGap: "
	    //          << trialElasticCircumGap << endl;
	    //     Info << "torsionCircumGap: "
	    //          << torsionCircumGap << endl;
	    //     Info << "twistCircumGap: "
	    //          << twistCircumGap << endl;
	    //     Info << "maxFrictionForce: "
	    //          << maxFrictionForce << endl;
	    //     Info << "forceCorrection: "
	    //          << forceCorrection << endl;
	    //     Info << "trialCircumContactForce: "
	    //          << mag(trialCircumContactForce) << endl;
	    //     Info << "normalContactForce: "
	    //          << mag(normalContactForce_) << endl;
	    // }

            trialElasticCircumGap = 0; // derivative should be 0
        }
        else if (mag(trialCircumContactForce) < SMALL)
        {
	    if (mag(trialElasticCircumGap)<SMALL)
	    {
                circumferentialContactForce_ = vector::zero;
	    }
	    else
	    {
		vector newCircumContactForce =
		    maxFrictionForce*circumDir
		   *sign(trialElasticCircumGap);

		circumferentialContactForce_ +=
		    (newCircumContactForce - circumferentialContactForce_);
	    }
	    // Info << "x3" << endl;
        }
        else
        {
            circumferentialContactForce_ = trialCircumContactForce;
            circumferentialSlipGap_ = oldCircumferentialSlipGap_;

	    // if
	    // (
	    //     (bI == 0)
	    //  && (segI == debugSeg)
	    //  && debug
            // )
	    // {
	    //     Info << "x2" << endl;
	    //     Info << circumferentialGap_ << endl;
	    //     Info << circumferentialSlipGap_ << endl;
	    //     Info << oldCircumferentialSlipGap_ << endl;
	    //     Info << "trialElasticCircumGap: "
	    //          << trialElasticCircumGap << endl;
	    //     Info << "torsionCircumGap: "
	    //          << torsionCircumGap << endl;
	    //     Info << "twistCircumGap: "
	    //          << twistCircumGap << endl;
	    //     Info << "maxFrictionForce: "
	    //          << maxFrictionForce << endl;
	    //     Info << "trialCircumContactForce: "
	    //          << mag(trialCircumContactForce) << endl;
	    //     Info << "normalContactForce: "
	    //          << mag(normalContactForce_) << endl;
	    // }
        }
	
	scalar rf = 1;
	circumferentialContactForce_ =
	    rf
	   *(
	        circumferentialContactForce_
	      - prevCircumferentialContactForce
	    )
	  + prevCircumferentialContactForce;

	circumferentialContactMoment_ = //vector::zero;
  	    sign(circumferentialContactForce_ & circumDir)*(-1)
	   *mag(circumferentialContactForce_)*beam.R(bI)*T; // ZT: check this (bI or nbI)

	// circumferentialContactForce_ = vector::zero;

	// Circumferential force derivative

        circumferentialContactForceDerivative_ = //tensor::zero;
  	    sign(twistCircumGap)*twistSgn*(-1)
           *frictionModel.contactForceDerivative
            (
                mag(trialElasticCircumGap)
            )
	   *((circumDir*T) & (spinTensor(n) & (I-(n*n))))
	   *(beam.R(bI)/dist); // ZT: check this (bI or nbI)

	// vector Tm = (T+neiT);
	// Tm /= mag(Tm) + SMALL;

	circumferentialContactForceThetaDerivative_ = //tensor::zero;
	    sign(torsionCircumGap)
	    // sign(trialElasticCircumGap)
	   // *sign(circumferentialContactForce_ & circumDir)
           *frictionModel.contactForceDerivative
            (
                mag(trialElasticCircumGap)
            )
	   *beam.R(bI)*(circumDir*T); // ZT: check this (bI or nbI)

	circumferentialContactMomentDerivative_ = //tensor::zero;
	    sign(torsionCircumGap)*(-1)
	    // sign(trialElasticCircumGap)
	   // *sign(circumferentialContactForce_ & circumDir)
           *frictionModel.contactForceDerivative
            (
                mag(trialElasticCircumGap)
            )
	   *sqr(beam.R(bI)) // ZT: check this (bI or nbI)
	   *(spinTensor(R) & (circumDir*T));


	// if
	// (
	//     (bI == 0)
	//  && (segI == debugSeg)
	//  && debug
	// )
	// {
	//     Info << "z2" << endl;
	//     Info << "twistCircumGap: " << twistCircumGap << endl;
	//     Info << "torsionCircumGap: " << torsionCircumGap << endl;
	//     Info << "circumferentialGap: "
	// 	 << circumferentialGap_ << endl;
	//     Info << "trialElasticCircumGap: "
	// 	 << trialElasticCircumGap << endl;
	//     Info << "sgn trialElasticCircumGap: "
	// 	 << sign(trialElasticCircumGap) << endl;
	//     Info << "trialCircumContactForce: "
	// 	 << trialCircumContactForce << endl;
	//     Info << "circumContactForce: "
	// 	 << circumferentialContactForce_ << endl;
	//     Info << "circumContactMoment: "
	// 	 << circumferentialContactMoment_ << endl;
	//     Info << "normalContactForce: "
	// 	 << mag(normalContactForce_) << endl;
	// }
	
        // // Tangential force derivative

        // vector dRdZeta = splines[bI].paramFirstDerivative(segI, zeta);
        // vector neiDRdZeta = splines[nbI].paramFirstDerivative(neiSegI, neiZeta);

        // vector dR2dZeta2 = splines[bI].paramSecondDerivative(segI, zeta);
        // vector neiDR2dZeta2 =
	//     splines[nbI].paramSecondDerivative(neiSegI, neiZeta);

        // vector delta = dist*n;

        // scalarSquareMatrix M(2, 0.0);
        // M[0][0] = (dRdZeta & neiDRdZeta);
        // M[0][1] = ((delta & neiDR2dZeta2) - (neiDRdZeta & neiDRdZeta));
        // M[1][0] = ((delta & dR2dZeta2) + (dRdZeta & dRdZeta));
        // M[1][1] = -(neiDRdZeta & dRdZeta);

        // scalarSquareMatrix invM = M.LUinvert();

        // vector Tm =
        //    -(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta)
        //    *sign(firstBeamGap_);

        // firstBeamTangContactForceDerivative_ =
        //     frictionModel.contactForceDerivative
        //     (
        //         mag(trialElasticTangGap)
        //     )*(T*Tm)*sign(trialElasticTangGap)*J;

        // vector neiTm =
        //    -(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta)
        //    *sign(secondBeamGap_);

        // secondBeamTangContactForceDerivative_ =
        //     frictionModel.contactForceDerivative
        //     (
        //         mag(neiTrialElasticTangGap)
        //     )*(neiT*neiTm)*sign(neiTrialElasticTangGap)*neiJ;
	
        // Info << bI << ", " <<  segI << ", "
	//      << nbI << ", " << neiSegI << ", " << neiZeta << ", "
	//      << mag(normalContactForce_) << ", "
	//      << normalGap_ << endl;
    }
    else
    {
        // Info << "else: " << bI << ", " <<  segI << ", "
	//      << nbI << ", " << neiSegI << ", " << neiZeta << ", "
	//      << mag(normalContactForce_) << ", "
	//      << normalGap_ << endl;

        axialContactForce_ = vector::zero;
        circumferentialContactForce_ = vector::zero;
        axialContactForceDerivative_ = tensor::zero;
        circumferentialContactForceDerivative_ = tensor::zero;
    }
}

//- SB: Returns weight value for  line contact in ABC contact
//   formulation (adopted from Meier's thesis, eqn. 4.62, pg 185)
Foam::scalar Foam::lineContact::allAngleLineContactWeight
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines,
    const scalar lowerContactAngleLimit,
    const scalar upperContactAngleLimit	  
) 
{
    scalar lineWeight = GREAT;

    label bI = firstBeam();
    label segI = firstBeamSegment();
    // scalar zeta = firstBeamZeta();  // not needed because contact is 
				       // evaluated at centre of cell

    label nbI = secondBeam();
    label neiSegI = secondBeamSegment();
    scalar neiZeta = secondBeamZeta();

    // scalar dist = GREAT;

    if (neiSegI > -1)
    {
	
	// Calculating contact angle for line contact
	// at zeta = 0 location (centre of CV)
        vector tan0 = splines[bI].firstDerivative(segI, 0);
        
	tan0 /= mag(tan0) + SMALL;

	vector tan1 =
		splines[nbI].firstDerivative(neiSegI, neiZeta);
			
	tan1 /= mag(tan1) + SMALL;

	scalar contactAngle = ::acos(mag(tan0 & tan1))*180.0/M_PI;
	
	//Info << "contact Angle " << contactAngle << endl;
	
	if (contactAngle < lowerContactAngleLimit)
	{
	    lineWeight = 1;
	}
	else if (contactAngle > upperContactAngleLimit)
	{
	    lineWeight = 0;
	}
	else
	{
	    scalar z = ::cos(contactAngle*M_PI/180.0);
	    scalar z1 = ::cos(lowerContactAngleLimit*M_PI/180.0);
	    scalar z2 = ::cos(upperContactAngleLimit*M_PI/180.0);
	    
	    scalar a = (M_PI*((z - z2)/(z1 - z2)));
	    
	    lineWeight = 0.5*(1 - ::cos(a));
	
	}
    }
    
    return lineWeight;
}



// This member function is probably not being used
Foam::scalar Foam::lineContact::calcNormalGap
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines
) const
{
    scalar curNormalGap = 0;

    label bI = firstBeam();
    label segI = firstBeamSegment();
    scalar zeta = firstBeamZeta();

    label nbI = secondBeam();
    label neiSegI = secondBeamSegment();
    scalar neiZeta = secondBeamZeta();

    // scalar dist = GREAT;

    if (neiSegI > -1)
    {
        vector curContactPoint = 
            splines[bI].position(segI, zeta);

        vector curNeiContactPoint = 
            splines[nbI].position(neiSegI, neiZeta);

        scalar dist = mag(curContactPoint - curNeiContactPoint);

        curNormalGap = dist - beam.R(bI) - beam.R(nbI);
    }
    
    return curNormalGap;
}


// ************************************************************************* //
