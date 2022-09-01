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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointContact::pointContact()
:
    firstBeam_(-1),
    firstBeamSegment_(-1),
    initialFirstBeamSegment_(-1),
    firstBeamOldSegment_(-1),
    firstBeamZeta_(0),
    initialFirstBeamZeta_(0),
    firstBeamOldZeta_(0),
    firstBeamLowerSegment_(-1),
    firstBeamUpperSegment_(-1),
    firstBeamWeight_(0),
    secondBeam_(-1),
    secondBeamSegment_(-1),
    initialSecondBeamSegment_(-1),
    secondBeamOldSegment_(-1),
    secondBeamZeta_(0),
    initialSecondBeamZeta_(0),
    secondBeamOldZeta_(0),
    secondBeamLowerSegment_(-1),
    secondBeamUpperSegment_(-1),
    secondBeamWeight_(0),
    delta_(GREAT),
    normalContactForce_(vector::zero),
    normalContactForceOffset_(0),
    normalContactForceDerivative_(tensor::zero),
    firstBeamTangContactForce_(vector::zero),
    firstBeamTangContactForceDerivative_(tensor::zero),
    secondBeamTangContactForce_(vector::zero),
    secondBeamTangContactForceDerivative_(tensor::zero),
    firstBeamGap_(0),
    firstBeamOldGap_(0),
    secondBeamGap_(0),
    secondBeamOldGap_(0),
    firstBeamSlipGap_(0),
    firstBeamOldSlipGap_(0),
    secondBeamSlipGap_(0),
    secondBeamOldSlipGap_(0),
    curTimeIndex_(-1)
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
    firstBeamZeta_(fbZeta),
    initialFirstBeamZeta_(0),
    firstBeamOldZeta_(0),
    firstBeamLowerSegment_(-1),
    firstBeamUpperSegment_(-1),
    firstBeamWeight_(0),
    secondBeam_(sb),
    secondBeamSegment_(sbSeg),
    initialSecondBeamSegment_(-1),
    secondBeamOldSegment_(-1),
    secondBeamZeta_(sbZeta),
    initialSecondBeamZeta_(0),
    secondBeamOldZeta_(0),
    secondBeamLowerSegment_(-1),
    secondBeamUpperSegment_(-1),
    secondBeamWeight_(0),
    delta_(GREAT),
    normalContactForce_(vector::zero),
    normalContactForceOffset_(0),
    normalContactForceDerivative_(tensor::zero),
    firstBeamTangContactForce_(vector::zero),
    firstBeamTangContactForceDerivative_(tensor::zero),
    secondBeamTangContactForce_(vector::zero),
    secondBeamTangContactForceDerivative_(tensor::zero),    
    firstBeamGap_(0),
    firstBeamOldGap_(0),
    secondBeamGap_(0),
    secondBeamOldGap_(0),
    firstBeamSlipGap_(0),
    firstBeamOldSlipGap_(0),
    secondBeamSlipGap_(0),
    secondBeamOldSlipGap_(0),
    curTimeIndex_(timeIndex)
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
    Info << "Inside pointContact::set() pointContact.C file" << endl;
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

        // Info << firstBeamOldSlipGap_ << ", "
        //      <<  secondBeamOldSlipGap_ << endl;
    }
    
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
    const bool agumentedLagrangian
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
	Info << "\nzeta " << zeta << " neiZeta " << neiZeta << endl;
	Info << "\nsegI " << segI << " neiSegI " << neiSegI << endl;
	
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
	
	Info << "fbLwSeg " << firstBeamLowerSegment_
	     << " fbUpSeg " << firstBeamUpperSegment_ << endl; 
	     
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
		if (zeta > 0.998)
		{
		    firstBeamWeight_ = 0.5;
		}
		else
		{
		    firstBeamWeight_ = (l-lx)/l;
		}
            }
            else
            {
		if (zeta < -0.998)
		{
		    firstBeamWeight_ = 0.5;
		}
		else
		{
		    firstBeamWeight_ = 1.0 - (l-lx)/l;
		}
            }            
        }
	Info << "firstBeamWeight " << firstBeamWeight_ << endl;
    }
    
    // Determine upper and lower segments for neighbour beam
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
		    
	Info << "sbLwSeg "  << secondBeamLowerSegment_
	     << " sbUpSeg " << secondBeamUpperSegment_ 
	     << endl;

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
		if (neiZeta > 0.998)
		{
		    secondBeamWeight_ =  0.5;
		}
		else
		{
		    secondBeamWeight_ = (l-lx)/l;
		}
            }
            else
            {
		if (neiZeta < -0.998)
		{
		    secondBeamWeight_ = 0.5;
		}
		else
		{
		     secondBeamWeight_ = 1.0 - (l-lx)/l;
		}
            }            
        }
	Info << "secondBeamWeight " << secondBeamWeight_ << endl;
    }

    // Calculate contact forces
    vector curContactPoint = 
        splines[bI].position(segI, zeta);

    vector curNeiContactPoint = 
        splines[nbI].position(neiSegI, neiZeta);

    scalar dist = mag(curContactPoint - curNeiContactPoint);

    vector n = curContactPoint - curNeiContactPoint;
    n /= mag(n) + SMALL;

    scalar gap = dist - beam.R(bI) - beam.R(nbI);

    delta() = dist;
    
    // SB: new pointWeight according to ABC formulation added
    scalar pointWeight =
	allAnglePointContactWeight
	(
	    beam,
	    splines,
	    lowerContactAngleLimit,
	    upperContactAngleLimit
	);
	
    // Info << "pointWeight " << pointWeight << endl;

    normalContactForce() = pointWeight*
        normalModel.contactForce
        (
            gap,
            normalContactForceOffset()
        )*n;
    // SB Added these two lines to print the normal cntact force value
    normalContactForce_ = pointWeight*
		normalModel.contactForce
		    (
			gap,
			normalContactForceOffset()
		    )*n;
		    
    Info << "\nNormal contact force " << pointWeight*normalContactForce_ 
	 << endl;

    normalContactForceDerivative() = pointWeight*
        normalModel.contactForceDerivative
        (
            gap,
            normalContactForceOffset()
        )*sqr(n);
    
    // Info << "sqr(n): " << sqr(n) << endl ;
    
    if (agumentedLagrangian)
    {
        if ((gap < -1e-4*(beam.R(bI) + beam.R(nbI))/2) || (gap > 0))
        {
            normalContactForceOffset() +=
            (    
                mag(normalContactForce())
              - normalContactForceOffset()
            );
        }
    }

    // Tangential force
    
    if
    (
        (mag(normalContactForce_) > SMALL)
     && (initialFirstBeamSegment_ == -1)
     && (firstBeamOldSegment_ == -1)
    )
    {
	 //   Info << "Inside tangential force calculation \n"
	  //   << "Line 360 pointContact.C file " << endl;
	     
        initialFirstBeamSegment_ = firstBeamSegment();
        initialFirstBeamZeta_ = firstBeamZeta();

        initialSecondBeamSegment_ = secondBeamSegment();
        initialSecondBeamZeta_ = secondBeamZeta();
	
        firstBeamOldSegment_ = firstBeamSegment();
        firstBeamOldZeta_ = firstBeamZeta();

        secondBeamOldSegment_ = secondBeamSegment();
        secondBeamOldZeta_ = secondBeamZeta();
    }

    if (mag(normalContactForce_) > SMALL)
    {
        scalar maxFrictionForce =
            frictionModel.maxFrictionForce
            (
                mag(normalContactForce())
            );
	
	// Info << "max friction force " << tab << maxFrictionForce << endl;
	
        label oldSegI = firstBeamOldSegment_;
        scalar oldZeta = firstBeamOldZeta_;
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
	firstBeamGap_ += firstBeamOldGap_;
	
        // label initSegI = initialFirstBeamSegment_;
        // scalar initZeta = initialFirstBeamZeta_;
        // scalar tangGap =
        //     splines[bI].length
        //     (
        //         segI,
        //         zeta
        //     )
        //   - splines[bI].length
        //     (
        //         initSegI,
        //         initZeta
        //     );

        vector T = splines[bI].paramFirstDerivative(segI, zeta);
        scalar J = mag(T);
        T /= J;

        scalar trialElasticTangGap = firstBeamGap_ - firstBeamOldSlipGap_;
        // scalar trialElasticTangGap = tangGap - firstBeamOldSlipGap_;

        vector trialTangContactForce = 
            frictionModel.contactForce
            (
                mag(trialElasticTangGap)
            )*T*sign(trialElasticTangGap);

        // Info << "g: " << tangGap << endl;

        if (mag(trialTangContactForce) > maxFrictionForce)
        { // Info << "line 430" << endl;
            firstBeamTangContactForce_ =
                maxFrictionForce*T*sign(trialElasticTangGap);

            firstBeamSlipGap_ =
                firstBeamOldSlipGap_
              + frictionModel.slipGapCorrection
                (
                    mag(trialElasticTangGap),
                    (sign(trialElasticTangGap)*T)
                  & (trialTangContactForce - firstBeamTangContactForce_)
                )*sign(trialElasticTangGap);

            trialElasticTangGap = firstBeamGap_ - firstBeamSlipGap_;
            // trialElasticTangGap = tangGap - firstBeamSlipGap_;
            trialElasticTangGap = 0; // Zero contact force derivative            
            // Info << "first beam combined" << endl;
        }
        else if (mag(trialTangContactForce) < SMALL)
        {
            firstBeamTangContactForce_ =
                maxFrictionForce*T*sign(trialElasticTangGap);
             
            firstBeamSlipGap_ = firstBeamGap_;
            trialElasticTangGap = firstBeamGap_ - firstBeamSlipGap_;
            // firstBeamSlipGap_ = tangGap;
            // trialElasticTangGap = tangGap - firstBeamSlipGap_;
        }
        else
        {
            firstBeamTangContactForce_ = trialTangContactForce;

            firstBeamSlipGap_ = firstBeamOldSlipGap_;
        }

        // Info << "gS: " << firstBeamSlipGap_ << endl;
        // Info << "gE: " << trialElasticTangGap << endl;
        // Info << "Ft: " << firstBeamTangContactForce_ << endl;
        
        // Info << firstBeamSlipGap_ << endl;
        
        label neiOldSegI = secondBeamOldSegment_;
        scalar neiOldZeta = secondBeamOldZeta_;
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
	secondBeamGap_ += secondBeamOldGap_;
	
        // label initNeiSegI = initialSecondBeamSegment_;
        // scalar initNeiZeta = initialSecondBeamZeta_;
        // scalar neiTangGap =
        //     splines[nbI].length
        //     (
        //         neiSegI,
        //         neiZeta
        //     )
        //   - splines[nbI].length
        //     (
        //         initNeiSegI,
        //         initNeiZeta
        //     );

        vector neiT = splines[nbI].paramFirstDerivative(neiSegI, neiZeta);
        scalar neiJ = mag(neiT);
        neiT /= neiJ;

        scalar neiTrialElasticTangGap = secondBeamGap_ - secondBeamOldSlipGap_;
        // scalar neiTrialElasticTangGap = neiTangGap - secondBeamOldSlipGap_;

        vector neiTrialTangContactForce =
            frictionModel.contactForce
            (
                mag(neiTrialElasticTangGap)
            )*neiT*sign(neiTrialElasticTangGap);

        if (mag(neiTrialTangContactForce) > maxFrictionForce)
        { //Info << "Line 514 " << endl;
            secondBeamTangContactForce_ =
                maxFrictionForce*neiT*sign(neiTrialElasticTangGap);

            secondBeamSlipGap_ =
                secondBeamOldSlipGap_
              + frictionModel.slipGapCorrection
                (
                    mag(neiTrialElasticTangGap),
                    (sign(neiTrialElasticTangGap)*neiT)
                  & (neiTrialTangContactForce - secondBeamTangContactForce_)
                )*sign(neiTrialElasticTangGap);          

            neiTrialElasticTangGap = secondBeamGap_ - secondBeamSlipGap_;
            // neiTrialElasticTangGap = neiTangGap - secondBeamSlipGap_;
            neiTrialElasticTangGap = 0; // derivative should be 0
            // Info << "second beam combined" << endl;
        }
        else if (mag(neiTrialTangContactForce) < SMALL)
        {
            secondBeamTangContactForce_ =
                maxFrictionForce*T*sign(neiTrialElasticTangGap);
             
            secondBeamSlipGap_ = secondBeamGap_;
            neiTrialElasticTangGap = secondBeamGap_ - secondBeamSlipGap_;
            // secondBeamSlipGap_ = neiTangGap;
            // neiTrialElasticTangGap = neiTangGap - secondBeamSlipGap_;
        }
        else
	{
            secondBeamTangContactForce_ = neiTrialTangContactForce;

            secondBeamSlipGap_ = secondBeamOldSlipGap_;
        }

        // Tangential force derivative
	// Info << "line 550 tangential force derivative " << endl;
        vector dRdZeta = splines[bI].paramFirstDerivative(segI, zeta);
        vector neiDRdZeta = splines[nbI].paramFirstDerivative(neiSegI, neiZeta);

        vector dR2dZeta2 = splines[bI].paramSecondDerivative(segI, zeta);
        vector neiDR2dZeta2 = splines[nbI].paramSecondDerivative(neiSegI, neiZeta);

        vector delta = dist*n;

        scalarSquareMatrix M(2, 0.0);
        M[0][0] = (dRdZeta & neiDRdZeta);
        M[0][1] = ((delta & neiDR2dZeta2) - (neiDRdZeta & neiDRdZeta));
        M[1][0] = ((delta & dR2dZeta2) + (dRdZeta & dRdZeta));
        M[1][1] = -(neiDRdZeta & dRdZeta);

        scalarSquareMatrix invM = M.LUinvert();

        vector Tm =
           -(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta)
           *sign(firstBeamGap_);
        // vector Tm =
        //    -(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta)
        //    *sign(tangGap);

        firstBeamTangContactForceDerivative_ =
            frictionModel.contactForceDerivative
            (
                mag(trialElasticTangGap)
            )*(T*Tm)*sign(trialElasticTangGap)*J;
      // + frictionModel.frictionCoeff()
      //  *normalModel.contactForceDerivative
      //   (
      //       gap,
      //       normalContactForceOffset()
      //   )*(T*n)*sign(trialElasticTangGap);

        vector neiTm =
           -(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta)
           *sign(secondBeamGap_);
        // vector neiTm =
        //    -(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta)
        //    *sign(neiTangGap);

        secondBeamTangContactForceDerivative_ =
            frictionModel.contactForceDerivative
            (
                mag(neiTrialElasticTangGap)
            )*(neiT*neiTm)*sign(neiTrialElasticTangGap)*neiJ;
      // - frictionModel.frictionCoeff()
      //  *normalModel.contactForceDerivative
      //   (
      //       gap,
      //       normalContactForceOffset()
      //   )*(neiT*n)*sign(trialElasticTangGap);
    }
    else
    {
	// Info << " Friction forces neglected \n"
	//     << " Tangential contact forces/derivatives set to zero \n"
	//     << " line 607 pointContact.C " << endl; 
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
