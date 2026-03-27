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

#include "pointContact.H"
#include "beamModel.H"
#include "scalarMatrices.H"
#include "beamContactModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointContact::pointContact()
:
    firstBeam_(-1),
    firstBeamSegment_(-1),
    initialFirstBeamSegment_(-1),
    firstBeamOldSegment_(-1),
    firstBeamPrevSegment_(-1),
    firstBeamZeta_(0),
    initialFirstBeamZeta_(0),
    firstBeamOldZeta_(0),
    firstBeamPrevZeta_(-2),
    firstBeamLowerSegment_(-1),
    firstBeamUpperSegment_(-1),
    firstBeamWeight_(0),
    secondBeam_(-1),
    secondBeamSegment_(-1),
    initialSecondBeamSegment_(-1),
    secondBeamOldSegment_(-1),
    secondBeamPrevSegment_(-1),
    secondBeamZeta_(0),
    initialSecondBeamZeta_(0),
    secondBeamOldZeta_(0),
    secondBeamPrevZeta_(-2),
    secondBeamLowerSegment_(-1),
    secondBeamUpperSegment_(-1),
    secondBeamWeight_(0),
    delta_(GREAT),
    normalGap_(GREAT),
    normalContactForce_(vector::zero),
    normalContactForceOffset_(0),
    normalContactForceDerivative_(tensor::zero),
    normalContactForceDn_(tensor::zero),
    firstBeamTangContactForce_(vector::zero),
    firstBeamTangContactForceDerivative_(tensor::zero),
    firstBeamTangContactForceDt_(tensor::zero),
    secondBeamTangContactForce_(vector::zero),
    secondBeamTangContactForceDerivative_(tensor::zero),
    secondBeamTangContactForceDt_(tensor::zero),
    couplingFirstBeamTangContactForceDerivative_(tensor::zero),
    couplingSecondBeamTangContactForceDerivative_(tensor::zero),
    couplingFirstBeamTangContactForceDt_(tensor::zero),
    couplingSecondBeamTangContactForceDt_(tensor::zero),
    firstBeamGap_(0),
    firstBeamOldGap_(0),
    firstBeamPrevGap_(0),
    secondBeamGap_(0),
    secondBeamOldGap_(0),
    secondBeamPrevGap_(0),
    firstBeamElasticGap_(0),
    secondBeamElasticGap_(0),
    firstBeamSlipGap_(0),
    firstBeamOldSlipGap_(0),
    firstBeamPrevSlipGap_(0),
    secondBeamSlipGap_(0),
    secondBeamOldSlipGap_(0),
    secondBeamPrevSlipGap_(0),
    curTimeIndex_(-1),
    alpha_(-1)
{}

Foam::pointContact::pointContact
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
    initialFirstBeamSegment_(-1),
    firstBeamOldSegment_(-1),
    firstBeamPrevSegment_(-1),
    firstBeamZeta_(fbZeta),
    initialFirstBeamZeta_(0),
    firstBeamOldZeta_(0),
    firstBeamPrevZeta_(-2),
    firstBeamLowerSegment_(-1),
    firstBeamUpperSegment_(-1),
    firstBeamWeight_(0),
    secondBeam_(sb),
    secondBeamSegment_(sbSeg),
    initialSecondBeamSegment_(-1),
    secondBeamOldSegment_(-1),
    secondBeamPrevSegment_(-1),
    secondBeamZeta_(sbZeta),
    initialSecondBeamZeta_(0),
    secondBeamOldZeta_(0),
    secondBeamPrevZeta_(-2),
    secondBeamLowerSegment_(-1),
    secondBeamUpperSegment_(-1),
    secondBeamWeight_(0),
    delta_(GREAT),
    normalGap_(GREAT),
    normalContactForce_(vector::zero),
    normalContactForceOffset_(0),
    normalContactForceDerivative_(tensor::zero),
    normalContactForceDn_(tensor::zero),
    firstBeamTangContactForce_(vector::zero),
    firstBeamTangContactForceDerivative_(tensor::zero),
    firstBeamTangContactForceDt_(tensor::zero),
    secondBeamTangContactForce_(vector::zero),
    secondBeamTangContactForceDerivative_(tensor::zero),
    secondBeamTangContactForceDt_(tensor::zero),
    couplingFirstBeamTangContactForceDerivative_(tensor::zero),
    couplingSecondBeamTangContactForceDerivative_(tensor::zero),
    couplingFirstBeamTangContactForceDt_(tensor::zero),
    couplingSecondBeamTangContactForceDt_(tensor::zero),
    firstBeamGap_(0),
    firstBeamOldGap_(0),
    firstBeamPrevGap_(0),
    secondBeamGap_(0),
    secondBeamOldGap_(0),
    secondBeamPrevGap_(0),
    firstBeamElasticGap_(0),
    secondBeamElasticGap_(0),
    firstBeamSlipGap_(0),
    firstBeamOldSlipGap_(0),
    firstBeamPrevSlipGap_(0),
    secondBeamSlipGap_(0),
    secondBeamOldSlipGap_(0),
    secondBeamPrevSlipGap_(0),
    curTimeIndex_(timeIndex),
    alpha_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointContact::set
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
    if (curTimeIndex_ < timeIndex)
    {
        firstBeamOldSegment_ = firstBeamSegment_;
        firstBeamOldZeta_ = firstBeamZeta_;

        secondBeamOldSegment_ = secondBeamSegment_;
        secondBeamOldZeta_ = secondBeamZeta_;

        firstBeamOldGap_ = firstBeamGap_; 
        secondBeamOldGap_ = secondBeamGap_;

        firstBeamOldSlipGap_ = firstBeamSlipGap_; 
        secondBeamOldSlipGap_ = secondBeamSlipGap_;
	
        curTimeIndex_  = timeIndex;
    }
    
    firstBeamPrevSegment_ = firstBeamSegment_;
    secondBeamPrevSegment_ = secondBeamSegment_;
    
    firstBeamPrevZeta_= firstBeamZeta_;
    secondBeamPrevZeta_ = secondBeamZeta_;
    
    firstBeamPrevGap_ = firstBeamGap_;
    secondBeamPrevGap_ = secondBeamGap_;
    
    firstBeamPrevSlipGap_= firstBeamSlipGap_;
    secondBeamPrevSlipGap_ = secondBeamSlipGap_;
    
    Info << "\nfbseg P " << firstBeamPrevSegment_
	 << " fbzeta P " << firstBeamPrevZeta_ 
	 << endl;
    
    firstBeam_ = fb;
    firstBeamSegment_ = fbSeg;
    firstBeamZeta_ = fbZeta;
    
    secondBeam_ = sb;
    secondBeamSegment_ = sbSeg;
    secondBeamZeta_ = sbZeta;
}


void Foam::pointContact::updateForce
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines,
    const normalContactModel& normalModel,
    const frictionContactModel& frictionModel,
    const scalar lowerContactAngleLimit,
    const scalar upperContactAngleLimit,
    const bool augmentedLagrangian
)
{
    Info << "\nInside pointContact::UpdateForce() in pointContact.C " << endl;
    
    label bI = firstBeam();
    label segI = firstBeamSegment();
    scalar zeta = firstBeamZeta();

    label nbI = secondBeam();
    label neiSegI = secondBeamSegment();
    scalar neiZeta = secondBeamZeta();
    

    // Determine upper and lower segments for first beam
    if (segI != -1)
    {	
        const label nSegments =
            beam.contact().splines()[bI].nSegments();
        
        if
        (
            (segI == 0 && zeta <= 0)
         || (segI == (nSegments-1) && zeta >= 0)
        )
        {
            firstBeamLowerSegment_ = segI;
            firstBeamUpperSegment_ = segI;
        }
        else if (zeta > 0)
        {
            firstBeamLowerSegment_ = segI;
            firstBeamUpperSegment_ = segI+1;
        }
        else
        {
            firstBeamLowerSegment_ = segI-1;
            firstBeamUpperSegment_ = segI;
        }
	
	     
	// Setting up the weights for first beam
        if (firstBeamLowerSegment_ == firstBeamUpperSegment_)
        {
            firstBeamWeight_ = 1.0;
        }
        else
        {
            scalar l =
                0.5*beam.contact().splines()[bI]
               .segLength(firstBeamLowerSegment_)
              + 0.5*beam.contact().splines()[bI]
               .segLength(firstBeamUpperSegment_);
	    
            scalar lx =
                mag(zeta)
               *beam.contact().splines()[bI]
               .segLength(segI)/2;

            if (segI == firstBeamLowerSegment_)
            {
		firstBeamWeight_ = (l-lx)/l;
            }
            else
            {
		firstBeamWeight_ = 1.0 - (l-lx)/l;
            }            
        }
    }
    
    // Determine upper and lower segments for the neighbour beam
    if (neiSegI != -1)
    {
        const label nSegments =
            beam.contact().splines()[nbI].nSegments();
        
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
    
    /*----------  PART I: NORMAL CONTACT FORCE CALCULATION   -----------*/
 
    // Calculate contact forces
    vector curContactPoint = 
        splines[bI].position(segI, zeta);

    vector curNeiContactPoint = 
        splines[nbI].position(neiSegI, neiZeta);
    
    Info << "segI: " << segI << "\nzeta: " << zeta
	 << "\nNeisegI: " << neiSegI << "\nnei zeta: " << neiZeta << endl;
    
    //   Info << "curContactPoint: " << curContactPoint
	 //   << "\nneiCurContactPoint: " << curNeiContactPoint << endl;
	 
    scalar dist = mag(curContactPoint - curNeiContactPoint);
    

    vector n = curContactPoint - curNeiContactPoint;
    n /= mag(n) + SMALL;

    scalar gap = dist - beam.R(bI) - beam.R(nbI);

    delta() = dist;
    
    normalGap_ = gap;
    

    Info << "gap: " << gap << endl;
    
    // SB: new pointWeight according to ABC formulation added
    scalar pointWeight =
	allAnglePointContactWeight
	(
	    beam,
	    splines,
	    lowerContactAngleLimit,
	    upperContactAngleLimit
	);
    
    Info << "normal contact force offset: " << normalContactForceOffset() << endl;

    // Explicit contact force contribution

	
    if (augmentedLagrangian)
    {
	//   //if ((gap < -1e-4*(beam.R(bI) + beam.R(nbI))/2) || (gap > 0))
	//   //if (gap < -1e-4*(beam.R(bI) + beam.R(nbI))/2)
	Info << "\nnormal contact offset " << normalContactForceOffset()
	    << endl;

	     
	if
	(
	    (
		(gap < -1e-3*(beam.R(bI) + beam.R(nbI))/2) || 
		(gap > 1e-3*(beam.R(bI) + beam.R(nbI))/2)
	    )
	     &&
	    (beam.iOuterCorr() > 1)
	)
	{
	    scalar lagrangeMultiplier = 0;
	    if 
	    (
		normalModel.epsilon()*gap < normalContactForceOffset()
	    )
	    {
		lagrangeMultiplier = normalModel.epsilon()*gap;
	    }
	    else
	    {
		lagrangeMultiplier = normalContactForceOffset();
	    }
	    
	    //normalContactForceOffset() +=
	    //(    
		//mag(normalContactForce())
	      //- normalContactForceOffset()
	    //);
	   
	    normalContactForceOffset() -= lagrangeMultiplier;
	}
	
	normalContactForce_ = pointWeight*
	(
	    normalContactForceOffset() +
	    normalModel.contactForce
	    (
		gap,
		0
	    )
	)*n; 

	Info << "\nnormal contact offset after update " << normalContactForceOffset()
	     << "\n" << endl;
    
    }
    else
    {

	normalContactForce_ = pointWeight*
	(
	    normalModel.contactForce
	    (
		gap,
		normalContactForceOffset_
	    )
	)*n;
    }
    
    Info << "normal contact force: " << normalContactForce_ << endl;
    
    // Implicit contact force derivative (Delta g) contribution
    normalContactForceDerivative_ = pointWeight*
	normalModel.contactForceDerivative
	(
	    gap,
	    normalContactForceOffset()
	)*sqr(n);
    
    // SB: Calculation of implicit contact force (Delta n) contribution
    tensor ImNN = tensor::I - n*n;
    
    vector dRdZeta = splines[bI].paramFirstDerivative(segI, zeta);
    vector neiDRdZeta = splines[nbI].paramFirstDerivative(neiSegI, neiZeta);
    vector dR2dZeta2 = splines[bI].paramSecondDerivative(segI, zeta);
    vector neiDR2dZeta2 = splines[nbI].paramSecondDerivative(neiSegI, neiZeta);

    scalarSquareMatrix M(2, 0.0);
	
    // Derivation available in point contact (frictionless) wriggers 1997 paper	     
    M[0][0] = (dRdZeta & neiDRdZeta);
    M[0][1] = (dist*(n & neiDR2dZeta2) - (neiDRdZeta & neiDRdZeta));
    M[1][0] = (dist*(n & dR2dZeta2) + (dRdZeta & dRdZeta));
    M[1][1] = -(neiDRdZeta & dRdZeta);

    // SB: LUinvert not available in ESI version
    // scalarSquareMatrix invM = M.LUinvert();

    // Calculating inverse by explicit logic
    scalar det = M[0][0]*M[1][1] - M[0][1]*M[1][0];
    if (mag(det) < SMALL)
        {
            FatalErrorInFunction
                << "Singular 2x2 contact matrix"
                << abort(FatalError);
        }
    
    scalarSquareMatrix invM(2, 0.0);
    invM[0][0] =  M[1][1]/det;
    invM[0][1] = -M[0][1]/det;
    invM[1][0] = -M[1][0]/det;
    invM[1][1] =  M[0][0]/det;
    
    tensor TT = 
	- (dRdZeta*(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta))
	+ (neiDRdZeta*(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta));
    
    // SB: Implicit contact force (Delta n) contribution
    normalContactForceDn_ = mag(normalContactForce_)
	*(
	    (ImNN + (ImNN & TT))/dist
	);
    
    //   if (
	    //   augmentedLagrangian && 
	    //   beam.iOuterCorr() > 1 
	    //   && 
	    //   (
		//   (gap < -1e-2*(beam.R(bI) + beam.R(nbI))/2) || 
		//   (gap > 1e-2*(beam.R(bI) + beam.R(nbI))/2)
	    //   )
	//   )
	//   {
	    //   scalar lagrangeMultiplier = 0;
	    
	    //   if 
	    //   (
		//   normalModel.epsilon()*gap < normalContactForceOffset()
	    //   )
	    //   {
		//   lagrangeMultiplier = normalModel.epsilon()*gap;
	    //   }
	    //   else
	    //   {
		//   lagrangeMultiplier = normalContactForceOffset();
	    //   }
	    
	    //   normalContactForceOffset_ -= lagrangeMultiplier;
	//   }
    
    
    /*----------  PART II: TANGENTIAL CONTACT FORCE (FRICTION) CALCULATION   -----------*/
    
    // Tangential force calculation
    
    /* In case the initialisation has not been done, ensure it here */
    if
    (
        (mag(normalContactForce_) > SMALL)
     && (initialFirstBeamSegment_ == -1)
     && (firstBeamOldSegment_ == -1)
    )
    {
        initialFirstBeamSegment_ = firstBeamSegment();
        initialFirstBeamZeta_ = firstBeamZeta();

        initialSecondBeamSegment_ = secondBeamSegment();
        initialSecondBeamZeta_ = secondBeamZeta();
	
        firstBeamOldSegment_ = firstBeamSegment();
        firstBeamOldZeta_ = firstBeamZeta();

        secondBeamOldSegment_ = secondBeamSegment();
        secondBeamOldZeta_ = secondBeamZeta();
    }
    
    //- SB added this 01st FEB 2023
    //   firstBeamOldGap_ = firstBeamGap_; 
    //   secondBeamOldGap_ = secondBeamGap_;

    //   firstBeamOldSlipGap_ = firstBeamSlipGap_; 
    //   secondBeamOldSlipGap_ = secondBeamSlipGap_;

    // Invoke friction only when beams in contact
    
    if (mag(normalContactForce_) > SMALL)
    {
	// Maximum friction allowed - mu*F_n
        scalar maxFrictionForce =
            frictionModel.maxFrictionForce
            (
                mag(normalContactForce())
            );
	//- Setting up the under-relaxation factor if outer iterations more than 10
	if (beam.iOuterCorr() > 8)
	{
	    alpha_ = beam.contact().alpha();
	    //   Info << "alpha : " << alpha_ << endl;
	}
	else
	{
	    alpha_ = 1;
	    //   Info << "alpha : " << alpha_ << endl;
	}
	
	// FIRST BEAM
        label oldSegI = firstBeamOldSegment_;
        scalar oldZeta = firstBeamOldZeta_;
	
	//   volVectorField W_prev = beam.solutionW().prevIter();
	//   volVectorField W = beam.solutionW();
	
	label globalSegI = beam.whichCell(bI, segI);
	
	//   /*----*/
	// Calculating the rotation matrix corresponding to old moving basis
	vector oldContactPoint = 
        splines[bI].position(oldSegI, oldZeta);
	
	vector oldNeiContactPoint = 
        splines[nbI].position(secondBeamOldSegment_, secondBeamOldZeta_);
	
	vector oldN = oldContactPoint - oldNeiContactPoint;
	oldN /= mag(oldN);
	
	vector oldT = splines[bI].paramFirstDerivative(oldSegI, oldZeta);
	
	vector oldB = oldT ^ oldN;
	
	vector e1{1,0,0};
	vector e2{0,1,0};
	vector e3{0,0,1};

	// Always initialise using curly brackets; do not assign using equal to
	// throws error
	tensor oldLambdaC
	    {
		oldT & e1, oldN & e1, oldB & e1,
		oldT & e2, oldN & e2, oldB & e2,
		oldT & e3, oldN & e3, oldB & e3
	    };
	    
	 // Calculating the change in the rotation tensor  direction
	 // due to change in normal of the contact point
	 vector T = splines[bI].paramFirstDerivative(segI, zeta);
	 scalar J = mag(T);
	 T /= J;
	 
	 vector B = T ^ n;
	 
	 // Increment rotation matrix for change in contact point
	 // This matrix is in spatial configuration
	 tensor delLambdaC
	    {
		T & oldT, n & oldT, B & oldT,
		T & oldN, n & oldN, B & oldN,
		T & oldB, n & oldB, B & oldB
	    };
	 
	 // Increment rotation matrix for change in contact point
	 // This matrix is in material configuration
	 tensor delLambda = oldLambdaC & (delLambdaC &  oldLambdaC.T());
	 
	 //   //Info << "del Lambda:\n" << delLambda << endl; 
	 
	 // This delLambda needs to be multiplied to the oldbeamGap
	 // most likely

	//   /*----*/
	
	// Previous iteration segI and zeta values
	label prevSegI = firstBeamPrevSegment_;
	scalar prevZeta = firstBeamPrevZeta_;
	
	// Calculate (scalar) increment in the first beam gap
	//   // For the first iteration, this increment is zero
        firstBeamGap_ =
            splines[bI].length
            (
                segI,
                zeta
            )
          - splines[bI].length
            (
                oldSegI,
                oldZeta
            );
    
	
	// Sign of the sliding tendency direction    
	scalar stilde = sign(firstBeamGap_);
	
	// Total gap increment vector in first beam
	vector fbGapIncrement =
	    splines[bI].position
	    (
		segI,
		zeta
	    )
	  - splines[bI].position
	   (
		oldSegI,
		oldZeta
	   );
	
	// (Scalar) tangential gap in the current configuration
	scalar fbtMag = mag(fbGapIncrement);
	firstBeamGap_ = firstBeamOldGap_ + stilde*fbtMag;
	//   firstBeamGap_ = firstBeamPrevGap_ + stilde*fbtMag;
	
	//   firstBeamGap_ += firstBeamOldGap_;

		
	// Direction in which the friction force should be applied
	vector fbT = fbGapIncrement/(fbtMag + SMALL);
	    
        //   // Gap increment vector in the current tangential direction
	//   fbGapIncrement_ -= (fbGapIncrement_ & n)*n;
    
	// Adding the increment to the old gap
	//   firstBeamGap_ += firstBeamOldGap_;	
	
	// Trial elastic tangential gap
        scalar trialElasticTangGap = firstBeamGap_ - firstBeamOldSlipGap_;
	 //   scalar trialElasticTangGap = firstBeamGap_ - firstBeamPrevSlipGap_;
	
        // scalar trialElasticTangGap = tangGap - firstBeamOldSlipGap_;

        //   vector trialTangContactForceVec = 
            //   frictionModel.contactForce
            //   (
                //   mag(trialElasticTangGap)
            //   )*fbT*stilde;
	    
	
	// Trial tangential friction force
	scalar trialTangContactForce = 
	    frictionModel.contactForce
	    (
		mag(trialElasticTangGap)
	    );
	
	if (false)
	{
	    Info << "First beam gap: " << firstBeamGap_
		 << "\nfbtMag: " << fbtMag << endl;
	    
	    Info  << "sign: " << stilde
		  << "\nT: " << T
		  << "\nfbT: " << fbT
		  << "\nfb trial gE: " << trialElasticTangGap << endl;
		
	    Info << "fb trial F:  " << trialTangContactForce << endl;
	}

	// Maximum tangential contact force can never exceed max friction force
	if (mag(trialTangContactForce) - maxFrictionForce  > SMALL)
	{	    
	    //   // (Explicit) tangential contact force
	    //   if (sTilde > 0)
	    //   {
		//   firstBeamTangContactForce_ = -stilde*maxFrictionForce*fbT;
	    //   }
	    //   else
	    //   {
		//   firstBeamTangContactForce_ = stilde*maxFrictionForce*fbT;
	    //   }
	    
	    // (Explicit) tangential contact force 
	    //   firstBeamTangContactForce_ = -stilde*maxFrictionForce*fbT; //actual expr
	    vector tmpFirstBeamForce = firstBeamTangContactForce_;
	    
	    firstBeamTangContactForce_ = 
		tmpFirstBeamForce + alpha_*(maxFrictionForce*fbT - tmpFirstBeamForce);
	    
	    //   firstBeamTangContactForce_ = maxFrictionForce*fbT;
	    
		 
	    // Slip gap correction as per litweka2000
	    scalar slipGapCorr = 
		frictionModel.slipGapCorrection
		(
		    mag(trialElasticTangGap), // Doubtful (not imp for 
		    // standard penalty method but have to check for regularised)
		    (
			mag(trialTangContactForce) 
			- mag(firstBeamTangContactForce_)
			//maxFrictionForce
			//   (trialTangContactForce - firstBeamTangContactForce_)&(fBt*sTilde)
		    )			    
		);
	    
	    // Updating the slip gap
	    firstBeamSlipGap_ = firstBeamOldSlipGap_ + slipGapCorr*stilde;
	    //   firstBeamSlipGap_ = firstBeamPrevSlipGap_ + slipGapCorr*stilde;
	    
	    // Updating the elastic gap  
            firstBeamElasticGap_ = firstBeamGap_ - firstBeamSlipGap_;
	    
	    //- Tangential force derivative (delta tangent gap) part
	    firstBeamTangContactForceDerivative_ =
		stilde*
		frictionModel.frictionCoeff()
	       *normalModel.contactForceDerivative
		(
		    gap,
		    normalContactForceOffset()
		)*(fbT*n);
	    
	    //- Tangential force (delta tangent) contribution
	    tensor B = - (dRdZeta*(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta));
	    
	    tensor ImTT = (tensor::I - fbT*fbT)/(fbtMag + SMALL);
	    
	    // THIS EXPR DOES NOT WORK
	    //   firstBeamTangContactForceDt_ = stilde*mag(firstBeamTangContactForce_)
		//   *(ImTT & (tensor::I + B));
	    
	    // THIS EXPR WORKS
	    firstBeamTangContactForceDt_ = 
		stilde*
		//   // maxFrictionForce*(ImTT & (B));
		mag(firstBeamTangContactForce_)*(ImTT & (B));
	    
	    couplingSecondBeamTangContactForceDerivative_ =
		firstBeamTangContactForceDerivative_;
		
	    couplingSecondBeamTangContactForceDt_ = 
		stilde*
		//   maxFrictionForce*(ImTT & B);
		mag(firstBeamTangContactForce_)*(ImTT & (B));
	    
	    
	    if (false)
	    {
		Info << "\nfirst beam status: slip " << endl;
		
		Info << "\nfb exp Force: " << firstBeamTangContactForce_ << endl;
		
		//   Info << "\nslipGapCorrection fb " << slipGapCorr 
		     //   << "\nfb old gS: " << firstBeamOldSlipGap_
		     //   << "\nfb new gS: " << firstBeamSlipGap_
		     //   << "\nfb gE: "  << trialElasticTangGap << endl;
		
		//   Info << "\nfb Force Der: " << firstBeamTangContactForceDerivative_
		     //   << "\nfb Force Dt: " << firstBeamTangContactForceDt_
		     //   << endl;
	    }
	    
	}
        else
        {
	    vector tmpFirstBeamForce = firstBeamTangContactForce_;
	    
	    firstBeamTangContactForce_ = 
		tmpFirstBeamForce + alpha_*(trialTangContactForce*fbT - tmpFirstBeamForce); 
		      
	    //   firstBeamTangContactForce_ = 
		//   trialTangContactForce*fbT;
	     //   -	(
		    //   frictionModel.contactForceDerivative
		    //   (
		    //   mag(trialElasticTangGap)
		    //   )*((fbT*fbT) & beam.solutionDW()[globalSegI])
		//   )
	     //   -	(
		    //   trialTangContactForce
		    //   *(
			//   (tensor::I - fbT*fbT)/(mag(fbGapIncrement) + SMALL)
		    //   )& beam.solutionDW()[globalSegI]
		//   );
	    
	    
	    // No increment in slip
            firstBeamSlipGap_ = firstBeamOldSlipGap_;
	    //   firstBeamSlipGap_ = firstBeamPrevSlipGap_;
	    
	    // Elastic gap
	    firstBeamElasticGap_ = firstBeamGap_ - firstBeamSlipGap_;

	    //- Tangential force derivative (delta tangent gap) part
	    tensor B = - (dRdZeta*(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta));
	    
	    //   vector Tm = - (invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta);
	    
	    tensor firstBeamTangContactDg = ((fbT * fbT) & (B));
	    
	    
	    // THIS EXPR DOES NOT WORK
	    //   tensor firstBeamTangContactDg = ((fbT * fbT) & (tensor::I + B)); 
	    
	    
	    //- Tangential force derivative (delta tangent gap)
	    firstBeamTangContactForceDerivative_ = 
		frictionModel.contactForceDerivative
		(
		    mag(firstBeamElasticGap_)
		)*firstBeamTangContactDg;
		// *(fbT*Tm);
	    
	    //- Tangential force (delta tangent) part
	    tensor ImTT = (tensor::I - fbT*fbT)/(fbtMag + SMALL);
	    
	    firstBeamTangContactForceDt_ = 
		////   trialTangContactForce*(ImTT & (B));
		mag(firstBeamTangContactForce_)*(ImTT & (B));
		
	    // THIS EXPR DOES NOT WORK
	    //- Tangential force (delta tangent) contribution
	    //   firstBeamTangContactForceDt_ = 
		//   mag(firstBeamTangContactForce_)
		//   *(ImTT & (tensor::I + B));
		

	    //- Coupling terms in second beam due to first beam
	    couplingSecondBeamTangContactForceDerivative_ = 
		(
		    frictionModel.contactForceDerivative
		    (
			mag(firstBeamElasticGap_)
		    )*((fbT * fbT) & B)
		    //   //*(fbT*Tm)
		);
		
	    couplingSecondBeamTangContactForceDt_ =
		(
		    //   trialTangContactForce*(ImTT & B)
		    mag(firstBeamTangContactForce_)*(ImTT & (B))
		);
	    
	    if (false)
	    {
		Info << "\nfirst beam status: stick " << endl;
		
		Info << "\nfb Force: " << firstBeamTangContactForce_ 
		     //   << "\nfb gE: " << trialElasticTangGap 
		     //   << "\n\nfb Force Der: " << firstBeamTangContactForceDerivative_
		     //   << "\nfb Force Dt " << firstBeamTangContactForceDt_
		     << endl;
	    }
        }
        
	
	// SECOND BEAM
        label neiOldSegI = secondBeamOldSegment_;
        scalar neiOldZeta = secondBeamOldZeta_;
	
	label globalNeiSegI = beam.whichCell(nbI, neiSegI);
     	    
	/*-----*/
	// Calculating the rotation matrix corresponding to old moving basis

	vector oldNeiT = splines[nbI].paramFirstDerivative(neiOldSegI, neiOldZeta);
	
	vector oldNeiB = oldNeiT ^ oldN;

	// Always initialise using curly brackets; do not assign using equal to
	// throws error
	tensor oldNeiLambdaC
	    {
		oldNeiT & e1, oldN & e1, oldNeiB & e1,
		oldNeiT & e2, oldN & e2, oldNeiB & e2,
		oldNeiT & e3, oldN & e3, oldNeiB & e3
	    };	    
	     
	 // Calculating the change in the rotation tensor  direction
	 // due to change in normal of the contact point
	 vector neiT = splines[nbI].paramFirstDerivative(neiSegI, neiZeta);
	 scalar neiJ = mag(neiT);
	 neiT /= neiJ;
	 
	 vector neiB = neiT ^ n;
	 
	 // Increment rotation matrix for change in contact point
	 // This matrix is in spatial configuration
	 tensor delNeiLambdaC
	    {
		neiT & oldNeiT, n & oldNeiT, neiB & oldNeiT,
		neiT & oldN,    n & oldN,    neiB & oldN,
		neiT & oldNeiB, n & oldNeiB, neiB & oldNeiB
	    };
	 
	 // Increment rotation matrix for change in contact point
	 // This matrix is in material configuration
	 tensor delNeiLambda = oldNeiLambdaC & (delNeiLambdaC &  oldNeiLambdaC.T());
	 
	 //   Info << "del nei Lambda:\n" << delNeiLambda << endl; 
	 
	 // This delNeiLambda needs to be multiplied to the oldbeamGap
	 // most likely
	 
	/*-----*/
	
	// Previous iteration segI and zeta values
	label neiPrevSegI = secondBeamPrevSegment_;
	scalar neiPrevZeta = secondBeamPrevZeta_;
	
	
       secondBeamGap_ =
	    splines[nbI].length
	    (
		neiSegI,
		neiZeta
	    )
	  - splines[nbI].length
	    (
		neiOldSegI,
		neiOldZeta
	    );
	
	// Sign of the sliding tendency direction    
	scalar neiStilde = sign(secondBeamGap_);
	
	// Gap increment vector in first beam
	vector sbGapIncrement =
	    splines[nbI].position
	    (
		neiSegI,
		neiZeta
	    )
	  - splines[nbI].position
	   (
		neiOldSegI,
		neiOldZeta
	   );
	
	// (Scalar) tangential gap in the current configuration
	scalar sbtMag = mag(sbGapIncrement);
	
	secondBeamGap_ = secondBeamOldGap_ + neiStilde*sbtMag;
	//   secondBeamGap_ = secondBeamPrevGap_ + neiStilde*sbtMag;
	//   secondBeamGap_ += secondBeamOldGap_;
	
	// Direction in which the friction force should be applied
	vector sbT = sbGapIncrement/(sbtMag + SMALL);
    
	// Trial elastic tangential gap
        scalar neiTrialElasticTangGap = secondBeamGap_ - secondBeamOldSlipGap_;
	//   scalar neiTrialElasticTangGap = secondBeamGap_ - secondBeamPrevSlipGap_;
	
	// Trial tangential friction force
	scalar neiTrialTangContactForce = 
	    frictionModel.contactForce
	    (
		mag(neiTrialElasticTangGap)
	    );
	
	//   vector neiTrialTangContactForceVec = 
            //   frictionModel.contactForce
            //   (
                //   mag(neiTrialElasticTangGap)
            //   )*sbT*neiStilde;
	
	if (false)
	{    
	    Info << "Second beam gap: " << secondBeamGap_ 
		 << "\nsbtMag: " << sbtMag << endl;
	    
	    Info << "sign: " << neiStilde
		 << "\nneiT: " << neiT
		 << "\nsbT: " << sbT
		 << "\nsb trial gE: " << neiTrialElasticTangGap << endl;
		
	    Info << "sb trial F:  " << neiTrialTangContactForce << endl;
	}

        if (mag(neiTrialTangContactForce) - maxFrictionForce > SMALL)
        { 
	    
	    //   if (neiStilde > 0)
	    //   {
		//   secondBeamTangContactForce_ = neiStilde*maxFrictionForce*sbT;
	    //   }
	    //   else
	    //   {
		//   secondBeamTangContactForce_ = -neiStilde*maxFrictionForce*sbT;
	    //   }
	    
	    vector tmpSecondBeamForce = secondBeamTangContactForce_;
	    
	    //   secondBeamTangContactForce_ = maxFrictionForce*sbT;
	    secondBeamTangContactForce_ =
		tmpSecondBeamForce + alpha_*(maxFrictionForce*sbT - tmpSecondBeamForce);
	    
	    // Slip gap correction as per litweka2000
	    scalar neiSlipGapCorr = 
		frictionModel.slipGapCorrection
		    (
			mag(neiTrialElasticTangGap),
			(
			    mag(neiTrialTangContactForce) 
			    //   - maxFrictionForce
			    - mag(secondBeamTangContactForce_)
			    //   (neiTrialTangContactForce 
			    //   - secondBeamTangContactForce_)&(sbT*neiStilde)
			)			    
		    );
				
	    // Updating the slip gap
	    secondBeamSlipGap_ = secondBeamOldSlipGap_ + neiStilde*neiSlipGapCorr;
	    //   secondBeamSlipGap_ = secondBeamPrevSlipGap_ + neiStilde*neiSlipGapCorr;
		 
	    // Updating the elastic gap
	    secondBeamElasticGap_ = secondBeamGap_ - secondBeamSlipGap_;
  
	    //- Tangential force derivative (delta tangent gap) part
	    secondBeamTangContactForceDerivative_ =
		neiStilde*
		    frictionModel.frictionCoeff()
		  *normalModel.contactForceDerivative
		    (
			gap,
			normalContactForceOffset()
		    )*(sbT*n);
	    
	    //- Tangential force (delta tangent) contribution
	    tensor C = (neiDRdZeta*(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta));
	    
	    tensor neiImTT = (tensor::I - sbT*sbT)/(sbtMag + SMALL);
	    
	    // THIS EXPR DOES NOT WORK
	    //   secondBeamTangContactForceDt_ = neiStilde*mag(secondBeamTangContactForce_)
		//   *(neiImTT & (tensor::I + C));
		
	    secondBeamTangContactForceDt_ = 
		neiStilde*
		////   maxFrictionForce*(neiImTT & (C));
		mag(secondBeamTangContactForce_)*(neiImTT & (C));
	    
	    //- Coupling terms for first beam due to second beam
	    couplingFirstBeamTangContactForceDerivative_ =
		secondBeamTangContactForceDerivative_;
	    
	    couplingSecondBeamTangContactForceDt_ =
		(
		    neiStilde*
		    //   maxFrictionForce*(neiImTT & C)
		    mag(secondBeamTangContactForce_)*(neiImTT & C)
		);
		
	    if (false)
	    {
		Info << "\nsecond beam status: slip " << endl;
		
		Info << "\nsb exp Force: " << secondBeamTangContactForce_ << endl;
		     // << "slipGapCorrection sb " << neiSlipGapCorr << endl;
		
		//   Info << "sb old gS: " << secondBeamOldSlipGap_
		     //   << "\nsb new gS: " << secondBeamSlipGap_ 
		     //   << "\nsb gE: "  << neiTrialElasticTangGap 
		     //   << "\n\nsb Force Der: " << secondBeamTangContactForceDerivative_
		     //   << "\nsb Force Dt: " << secondBeamTangContactForceDt_
		     //   << endl;
	    }

        }
        else
	{
	    vector tmpSecondBeamForce = secondBeamTangContactForce_;
	    
	    secondBeamTangContactForce_ =
		tmpSecondBeamForce + alpha_*(neiTrialTangContactForce*sbT - tmpSecondBeamForce);
			    
	    //   secondBeamTangContactForce_ = 
		//   neiTrialTangContactForce*sbT;
	     //   -	(
		    //   frictionModel.contactForceDerivative
		    //   (
		    //   mag(neiTrialElasticTangGap)
		    //   )*((sbT*sbT) & beam.solutionDW()[globalNeiSegI])
		//   )
	     //   -	(
		    //   neiTrialTangContactForce
		    //   *(
			//   (tensor::I - sbT*sbT)/(mag(sbGapIncrement) + SMALL)
		    //   )& beam.solutionDW()[globalNeiSegI]
		//   );
	    
	    
	    // No increment in slip gap
            secondBeamSlipGap_ = secondBeamOldSlipGap_;
	    //   secondBeamSlipGap_ = secondBeamPrevSlipGap_;
	    
	    // Elastic Gap
	    secondBeamElasticGap_ = secondBeamGap_ - secondBeamSlipGap_;
	    

	    //- Tangential force derivative (delta tangent gap) part
	    tensor C = (neiDRdZeta*(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta));
	    
	    //   vector neiTm = (invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta);
	    
	    // THIS EXPR DOES NOT WORK
	    //   tensor secondBeamTangContactDg = ((sbT * sbT) & (tensor::I + C));
	    
	    //- Tangential force (delta tangent gap) contribution
	    tensor secondBeamTangContactDg = ((sbT * sbT) & (C));
	    
	    secondBeamTangContactForceDerivative_ = 
		frictionModel.contactForceDerivative
		(
		    mag(secondBeamElasticGap_)
		)*secondBeamTangContactDg;
		//*(sbT*neiTm);
	    
	    
	    
	    //- Tangential force (delta tangent) contribution
	    tensor neiImTT = (tensor::I - sbT*sbT)/(sbtMag + SMALL);
	    
	    // THIS EXPR DOES NOT WORK
	    //   secondBeamTangContactForceDt_ = 
		//   mag(secondBeamTangContactForce_)
		//   *(neiImTT & (tensor::I + C));
		
	    // THIS WORKS 
	    secondBeamTangContactForceDt_ = 
		////   neiTrialTangContactForce*(neiImTT & (C));
		mag(secondBeamTangContactForce_)*(neiImTT & (C));


		
	    //- Coupling terms for first beam due to second beam 
	    couplingFirstBeamTangContactForceDerivative_ = 
		(
		    frictionModel.contactForceDerivative
		    (
			mag(secondBeamElasticGap_)
		    )*((sbT * sbT) & C)
		    //   //*(sbT*neiTm)
		    
		);
		
	    couplingFirstBeamTangContactForceDt_ = 
		(
		    //   neiTrialTangContactForce*(neiImTT & C)
		    mag(secondBeamTangContactForce_)*(neiImTT & (C))
		);
	    
	    if (false)
	    {
		Info << "\nsecond beam status: stick" << endl;
		Info << "\nsb exp Force: " << secondBeamTangContactForce_
		     //   << "\nsb gE: " << neiTrialElasticTangGap 
		     //   << "\n\nsb Force Der: " << secondBeamTangContactForceDerivative_
		     //   << "\nsb Force Dt: " << secondBeamTangContactForceDt_
		     << endl;
	    }
		    
        }

        //   // Tangential force derivative
        //   vector dRdZeta = splines[bI].paramFirstDerivative(segI, zeta);
        //   vector neiDRdZeta = splines[nbI].paramFirstDerivative(neiSegI, neiZeta);

        //   vector dR2dZeta2 = splines[bI].paramSecondDerivative(segI, zeta);
        //   vector neiDR2dZeta2 = splines[nbI].paramSecondDerivative(neiSegI, neiZeta);

        //   vector delta = dist*n;

        //   scalarSquareMatrix M(2, 0.0);
        //   M[0][0] = (dRdZeta & neiDRdZeta);
        //   M[0][1] = ((delta & neiDR2dZeta2) - (neiDRdZeta & neiDRdZeta));
        //   M[1][0] = ((delta & dR2dZeta2) + (dRdZeta & dRdZeta));
        //   M[1][1] = -(neiDRdZeta & dRdZeta);

        //   scalarSquareMatrix invM = M.LUinvert();

        //   vector Tm =
           //   -(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta)
           //   *sign(firstBeamGap_);
        // // vector Tm =
        // //    -(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta)
        // //    *sign(tangGap);

	/*
        firstBeamTangContactForceDerivative_ =
            frictionModel.contactForceDerivative
            (
                mag(trialElasticTangGap)
            )*(fBt*Tm)*sign(trialElasticTangGap);
        */
      //   // // + frictionModel.frictionCoeff()
      //   // //  *normalModel.contactForceDerivative
      //   // //   (
       //   // gap,
       //   // normalContactForceOffset()
      //   // //   )*(T*n)*sign(trialElasticTangGap);

        //   vector neiTm =
           //   -(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta)
           //   *sign(secondBeamGap_);
        //// vector neiTm =
        ////    -(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta)
        ////    *sign(neiTangGap);
	//   Info << "fb trial elastic gap "  << trialElasticTangGap << endl;
	//   Info << "sb trial elastic gap "  << neiTrialElasticTangGap << endl;

        //   secondBeamTangContactForceDerivative_ =
            //   frictionModel.contactForceDerivative
            //   (
                //   mag(neiTrialElasticTangGap)
            //   )*(sbTangConForceDir*neiTm)*sign(neiTrialElasticTangGap);
      //// - frictionModel.frictionCoeff()
      ////  *normalModel.contactForceDerivative
      ////   (
      ////       gap,
      ////       normalContactForceOffset()
      ////   )*(neiT*n)*sign(trialElasticTangGap);
      //   Info << "fb tang. con. force derivative " << firstBeamTangContactForceDerivative_ << endl;
      //   Info << "sb tang. con. force derivative " << secondBeamTangContactForceDerivative_ << endl;
    }
    else
    {
        firstBeamTangContactForce_ = vector::zero;
        secondBeamTangContactForce_ = vector::zero;
        firstBeamTangContactForceDerivative_ = tensor::zero;
        secondBeamTangContactForceDerivative_ = tensor::zero;
    }
}

//- SB: Returns weight value for point contact in ABC contact
//   formulation (adopted from Meier's thesis, eqn. 4.62, pg 185)
Foam::scalar Foam::pointContact::allAnglePointContactWeight
(
    const beamModel& beam,
    const PtrList<HermiteSpline>& splines,
    const scalar lowerContactAngleLimit,
    const scalar upperContactAngleLimit	  
) 
{
    scalar pointWeight = GREAT;

    label bI = firstBeam();
    label segI = firstBeamSegment();
    scalar zeta = firstBeamZeta();

    label nbI = secondBeam();
    label neiSegI = secondBeamSegment();
    scalar neiZeta = secondBeamZeta();

    // scalar dist = GREAT;

    if (neiSegI > -1)
    {
	
	// Calculating contact angle for point contact location
        vector tan0 = splines[bI].firstDerivative(segI, zeta);
        
	tan0 /= mag(tan0) + SMALL;

	vector tan1 =
		splines[nbI].firstDerivative(neiSegI, neiZeta);
			
	tan1 /= mag(tan1) + SMALL;

	scalar contactAngle = ::acos(mag(tan0 & tan1))*180.0/M_PI;
	
	// Info << "contact Angle " << contactAngle << endl;
	
	
	if (contactAngle < lowerContactAngleLimit)
	{
	    pointWeight = 0;
	}
	else if (contactAngle > upperContactAngleLimit)
	{
	    pointWeight = 1;
	}
	else
	{
	    scalar z = ::cos(contactAngle*M_PI/180.0);
	    scalar z1 = ::cos(lowerContactAngleLimit*M_PI/180.0);
	    scalar z2 = ::cos(upperContactAngleLimit*M_PI/180.0);
	    
	    scalar a = (M_PI*((z - z2)/(z1 - z2)));
	    
	    pointWeight = (
			    1.0 - 0.5*(1 - ::cos(a))
			  );
	
	}
    }
    
    return pointWeight;
}


// ************************************************************************* //
