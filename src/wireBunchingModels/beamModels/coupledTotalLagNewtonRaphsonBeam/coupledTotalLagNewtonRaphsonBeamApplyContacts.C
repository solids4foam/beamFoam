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
#include "scalarMatrices.H"
#include "spinTensor.H"
#include "beamHelperFunctions.H"
#include <Eigen/Sparse>
//#include "lagrangeMultipliers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

// SB (2026 - ESI):  New defn for off-diag coupling coeffs of beam contact
typedef Pair<scalarSquareMatrix> lcCoeffPair;
typedef List<lcCoeffPair> lcCoeffPairList;
typedef List<lcCoeffPairList> lcCoeffPairListList;

typedef Pair<scalarRectangularMatrix> pcCoeffPair;
typedef List<pcCoeffPair> pcCoeffPairList;
typedef List<pcCoeffPairList> pcCoeffPairListList;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void coupledTotalLagNewtonRaphsonBeam::applyLineContact
(
    Field<scalarSquareMatrix>& diag,
    Field<scalarSquareMatrix>& lower,
    Field<scalarSquareMatrix>& upper,
    Field<scalarRectangularMatrix>& source
    //    multibeamFvBlockMatrix& eqn
)
{
    if (mesh().cellZones().size() > 1)
    {
    	// Info << "Inside applyLineContact() in applyContacts.C \n" << endl;
        label nBeams = contact().splines().size();
        
        // // Get matrix diagonal
        // tensor6Field& diag = eqn.diag().asSquare();

        // // Grab source
        // vector6Field& source = eqn.source();

        // Grab beam coupling coeffs (foam-extend version)
        //        PtrList<tensor6PairListList>& lcCoeffs = eqn.lineContactCoeffs();

        // SB (2026 ESI) - A storage container to collect lcCoeffs for ESI version
        // IMPORTANT - However, these coefficients are simply collected but not actually
        // added to the block matrix of beam solver -
        // The lcCoeffs (and pcCoeffs) were passed to multiBeamFvMatrix class where
        // contributions were added. Since we do not use the class, the contribution of
        // the coeffs is simply being lost at this point.        
        PtrList<lcCoeffPairListList> lcCoeffs(nBeams);

        const lcCoeffPair lcZero
        (
            scalarSquareMatrix(3, Zero),
            scalarSquareMatrix(3, Zero)
        );

        for (label bI=0; bI<nBeams; ++bI)
        {
            const label nSeg = contact().splines()[bI].nSegments();

            lcCoeffs.set
            (
                bI,
                new lcCoeffPairListList(nSeg)
            );

            for (label segI=0; segI<nSeg; ++segI)
            {
                lcCoeffs[bI][segI].setSize(nBeams, lcZero);
            }
        }

        // Explicit part of the line contact force
        label start = 0;
        // Info << "Explicit line contact force: start " << endl;
        for (label bI=0; bI<nBeams; bI++)
        {
            const lineContactListList& curLineContacts =
                contact().lineContacts()[bI];

            vectorField curq
            (
                contact().splines()[bI].nSegments(),
                vector::zero
            );
            vectorField curm
            (
                contact().splines()[bI].nSegments(),
                vector::zero
            );

	    //-SB added
            const surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

	    const surfaceScalarField& dc = mesh().deltaCoeffs();
	    //-

            for
            (
                label segI=0;
                segI < contact().splines()[bI].nSegments();
                segI++
            )
            {
                label globalSegIndex = start + segI;

                label cellID =
                    localCellIndex(globalSegIndex);

                if (cellID != -1)
                {
                    for (label nbI=0; nbI<nBeams; nbI++)
                    {
                        if (nbI != bI) // No self-contact
                        {
			    //-SB added
			    const label neiSegI =
                                curLineContacts[segI][nbI].secondBeamSegment();

			    const scalar neiZeta =
				curLineContacts[segI][nbI].secondBeamZeta();

			    label globalNeiSegI = whichCell(nbI, neiSegI);
			    //-

                            const vector& curNormalContactForce =
                                curLineContacts[segI][nbI].normalContactForce();

                            if (mag(curNormalContactForce) > SMALL)
                            {
                                curq[segI] +=
                                    curLineContacts[segI][nbI]
                                   .normalContactForce()
                                  + curLineContacts[segI][nbI]
                                   .circumferentialContactForce()
                                  + curLineContacts[segI][nbI]
                                   .axialContactForce();

                                curm[segI] +=
                                    curLineContacts[segI][nbI]
                                   .circumferentialContactMoment();
                            }

			    //- SB added (16 Feb 2023)
			    vector DR = vector::zero;

			    if (neiZeta > 0)
			    {
                    //	label faceID = findIndex(mesh().owner(), globalNeiSegI);
                    label faceID = mesh().owner().find(globalNeiSegI);

				if (faceID == -1) // last cell
				{
				    //   Info << "bI " << bI << " nbI " << nbI
					//   << " neiSegI "<< neiSegI <<  " last cell DR " << endl;
				    // const unallocLabelList& faceCells =
                    const labelUList& faceCells =
					mesh().boundary()[endPatchIndex(nbI)].faceCells();

                    // label bFaceID = findIndex(faceCells, globalNeiSegI);
                    label bFaceID = faceCells.find(globalNeiSegI);

				    DR = neiZeta
				       *dRdS.boundaryField()[endPatchIndex(nbI)][bFaceID]
				       /dc.boundaryField()[endPatchIndex(nbI)][bFaceID];
				}
				else
				{
				    DR = 0.5*neiZeta*dRdS.internalField()[faceID]
				       /dc.internalField()[faceID];

				}
			    }
			    else
			    {
                    // label faceID = findIndex(mesh().neighbour(), globalNeiSegI);
                label faceID = mesh().neighbour().find(globalNeiSegI);
				if (faceID == -1) // first cell
				{
				    //   Info << "bI " << bI << " nbI " << nbI
					//   << " neiSegI "<< neiSegI <<  " first cell DR " << endl;
				    // const unallocLabelList& faceCells =
                    const labelUList& faceCells =
					mesh().boundary()[startPatchIndex(nbI)].faceCells();

                    // label bFaceID = findIndex(faceCells, globalNeiSegI);
                    label bFaceID = faceCells.find(globalNeiSegI);

				    DR = neiZeta
				       *dRdS.boundaryField()[startPatchIndex(nbI)][bFaceID]
				       /dc.boundaryField()[startPatchIndex(nbI)][bFaceID];
				}
				else
				{
				    DR = 0.5*neiZeta*dRdS.internalField()[faceID]
				       /dc.internalField()[faceID];

				}
			    }

			    //- Line Contact force acting at a point on the neighbouring
			    //- beam at neiZeta location
			    vector neiContactForce = - curNormalContactForce*L()[cellID];
			    vector neiContactMoment = (spinTensor(DR) & neiContactForce);

			    //   source[globalNeiSegI](0) =- neiContactForce.x();
			    //   source[globalNeiSegI](1) =- neiContactForce.y();
			    //   source[globalNeiSegI](2) =- neiContactForce.z();

			    //   source[globalNeiSegI](3) =- neiContactMoment.x();
			    //   source[globalNeiSegI](4) =- neiContactMoment.y();
			    //   source[globalNeiSegI](5) =- neiContactMoment.z();

			    //-
                        }
                    }
                }
            }

            forAll(curq, segI)
            {
                label globalSegIndex = start + segI;

                label cellID =
                    localCellIndex(globalSegIndex);

		//   if (mag(curq[segI]) != 0)
		//   {
		    //   if
		    //   (
		    //   segI == 0
		    //   || segI == (contact().splines()[bI].nSegments() - 1)
		    //   )
		    //   {
		    //   Info << "lc force " << mag(curq[segI]) << tab
			 //   << "segI " << segI << tab
			 //   << "cellID " << cellID << endl;
		    //   }
		//   }

                if (cellID != -1)
                {
                    // W equation
                    source[cellID](0, 0) -= curq[segI].x()*L()[cellID];
                    source[cellID](1, 0) -= curq[segI].y()*L()[cellID];
                    source[cellID](2, 0) -= curq[segI].z()*L()[cellID];

                    // Theta equation
                    source[cellID](3, 0) -= curm[segI].x()*L()[cellID];
                    source[cellID](4, 0) -= curm[segI].y()*L()[cellID];
                    source[cellID](5, 0) -= curm[segI].z()*L()[cellID];
                }
            }

            start += curq.size();
        }

        // Implicit part of the line contact force
        // Info << "Implicit line contact force: start" << endl;
        start = 0;
        for (label bI=0; bI<nBeams; bI++)
	{
            const lineContactListList& curLineContacts =
	    contact().lineContacts()[bI];

	    label nSeg = contact().splines()[bI].nSegments();

	    for (label segI = 0; segI < nSeg; segI++)
	    {
		label globalSegI = start + segI;

                label cellID = localCellIndex(globalSegI);

                if (cellID != -1)
                {
                    for (label nbI=0; nbI<nBeams; nbI++)
                    {
                        if (nbI != bI) // No self-contact
                        {
                            //   const label neiSegI =
                                //   curLineContacts[segI][nbI].secondBeamSegment();

                            const vector& curContactForce =
                                curLineContacts[segI][nbI].normalContactForce();

                            if (mag(curContactForce) > SMALL)
                            {
                                //   const scalar neiSegParam =
                                    //   curLineContacts[segI][nbI].secondBeamZeta();

                                // const label neiLowerSegI =
                                //     curLineContacts[segI][nbI]
                                //    .secondBeamLowerSegment();
                                // const label neiUpperSegI =
                                //     curLineContacts[segI][nbI]
                                //    .secondBeamUpperSegment();

                                //   const scalar& curContactDistance =
                                    //   curLineContacts[segI][nbI].delta();

                                tensor curContactForceDerivative =
                                    curLineContacts[segI][nbI]
                                   .normalContactForceDerivative()
                                   *L()[cellID];

                                //   const tensor curAxialContactForceDerivative =
                                    //   curLineContacts[segI][nbI]
                                   //   .axialContactForceDerivative();

                                //   const vector curCircumContactForce =
                                    //   curLineContacts[segI][nbI]
                                   //   .circumferentialContactForce();

                                //   const tensor curCircumContactForceDerivative =
                                    //   curLineContacts[segI][nbI]
                                   //   .circumferentialContactForceDerivative();

                                //   const tensor curCircumContactForceThetaDerivative =
                                    //   curLineContacts[segI][nbI]
                                   //   .circumferentialContactForceThetaDerivative();

                                //   const tensor curContactMomentDerivative =
                                    //   curLineContacts[segI][nbI]
                                   //   .circumferentialContactMomentDerivative();

                                //   const vector dRdZeta =
                                    //   contact().splines()[nbI]
                                   //   .paramFirstDerivative(neiSegI, neiSegParam);

                                // What is happening here? Check with ZT

                                //   tensor ImNN = tensor::zero;
                                //   vector N = curContactForce;
                                //   if (mag(N) > SMALL)
                                //   {
                                    //   N /= mag(N);
                                    //   ImNN = tensor::I - (N*N);
                                //   }

                                //   const tensor TT = //tensor::zero;
                                    //   (dRdZeta*dRdZeta)/((dRdZeta & dRdZeta));

                                //   curContactForceDerivative +=
                                //   (
				    //mag(curContactForce)*(ImNN - (ImNN & TT))
                                   //*L()[cellID]/curContactDistance

				    //   // mag(curContactForce)*(ImNN - (ImNN & TT))
                                   //   // /curContactDistance

				   //   //oiginal eqn is below
				   //   mag(curContactForce)*(ImNN - TT)
                                   //   *L()[cellID]/curContactDistance
                                //   );

                                //   // Add frictional component
                                //   if (mag(curCircumContactForce) > SMALL)
                                //   {
                                    //   curContactForceDerivative +=
                                        //   curCircumContactForceDerivative
                                       //   *L()[cellID];

                                    //   curContactForceDerivative +=
                                        //   curAxialContactForceDerivative
                                       //   *L()[cellID];
                                //   }

                                //- Contact force derivative (W)
                                diag[cellID](0,0) +=
                                    curContactForceDerivative.xx();
                                diag[cellID](0,1) +=
                                    curContactForceDerivative.xy();
                                diag[cellID](0,2) +=
                                    curContactForceDerivative.xz();

                                diag[cellID](1,0) +=
                                    curContactForceDerivative.yx();
                                diag[cellID](1,1) +=
                                    curContactForceDerivative.yy();
                                diag[cellID](1,2) +=
                                    curContactForceDerivative.yz();

                                diag[cellID](2,0) +=
                                    curContactForceDerivative.zx();
                                diag[cellID](2,1) +=
                                    curContactForceDerivative.zy();
                                diag[cellID](2,2) +=
                                    curContactForceDerivative.zz();

                                //   //- Contact force derivative (Theta)
                                //   diag[cellID](0,3) +=
                                    //   curCircumContactForceThetaDerivative.xx()
                                   //   *L()[cellID];
                                //   diag[cellID](0,4) +=
                                    //   curCircumContactForceThetaDerivative.xy()
                                   //   *L()[cellID];
                                //   diag[cellID](0,5) +=
                                    //   curCircumContactForceThetaDerivative.xz()
                                   //   *L()[cellID];

                                //   diag[cellID](1,3) +=
                                    //   curCircumContactForceThetaDerivative.yx()
                                   //   *L()[cellID];
                                //   diag[cellID](1,4) +=
                                    //   curCircumContactForceThetaDerivative.yy()
                                   //   *L()[cellID];
                                //   diag[cellID](1,5) +=
                                    //   curCircumContactForceThetaDerivative.yz()
                                   //   *L()[cellID];

                                //   diag[cellID](2,3) +=
                                    //   curCircumContactForceThetaDerivative.zx()
                                   //   *L()[cellID];
                                //   diag[cellID](2,4) +=
                                    //   curCircumContactForceThetaDerivative.zy()
                                   //   *L()[cellID];
                                //   diag[cellID](2,5) +=
                                    //   curCircumContactForceThetaDerivative.zz()
                                   //   *L()[cellID];

                                //   //- Contact moment derivative
                                //   diag[cellID](3,3) +=
                                    //   curContactMomentDerivative.xx()
                                   //   *L()[cellID];
                                //   diag[cellID](3,4) +=
                                    //   curContactMomentDerivative.xy()
                                   //   *L()[cellID];
                                //   diag[cellID](3,5) +=
                                    //   curContactMomentDerivative.xz()
                                   //   *L()[cellID];

                                //   diag[cellID](4,3) +=
                                    //   curContactMomentDerivative.yx()
                                   //   *L()[cellID];
                                //   diag[cellID](4,4) +=
                                    //   curContactMomentDerivative.yy()
                                   //   *L()[cellID];
                                //   diag[cellID](4,5) +=
                                    //   curContactMomentDerivative.yz()
                                   //   *L()[cellID];

                                //   diag[cellID](5,3) +=
                                    //   curContactMomentDerivative.zx()
                                   //   *L()[cellID];
                                //   diag[cellID](5,4) +=
                                    //   curContactMomentDerivative.zy()
                                   //   *L()[cellID];
                                //   diag[cellID](5,5) +=
                                    //   curContactMomentDerivative.zz()
                                   //   *L()[cellID];

                                // Off-diagonal

                                scalar w0 =
                                    curLineContacts[segI][nbI].secondBeamWeight();
                                scalar w1 = 1.0 - w0;

				//   Info << "weight w0: " << w0 << endl;
				//   Info << "w1: " << w1 << endl;

                                // Contact force derivative (W)

                                //- lower
                                lcCoeffs[bI][segI][nbI].first()(0, 0) -=
                                    curContactForceDerivative.xx()*w0;
                                lcCoeffs[bI][segI][nbI].first()(0, 1) -=
                                    curContactForceDerivative.xy()*w0;
                                lcCoeffs[bI][segI][nbI].first()(0, 2) -=
                                    curContactForceDerivative.xz()*w0;

                                lcCoeffs[bI][segI][nbI].first()(1, 0) -=
                                    curContactForceDerivative.yx()*w0;
                                lcCoeffs[bI][segI][nbI].first()(1, 1) -=
                                    curContactForceDerivative.yy()*w0;
                                lcCoeffs[bI][segI][nbI].first()(1, 2) -=
                                    curContactForceDerivative.yz()*w0;

                                lcCoeffs[bI][segI][nbI].first()(2, 0) -=
                                    curContactForceDerivative.zx()*w0;
                                lcCoeffs[bI][segI][nbI].first()(2, 1) -=
                                    curContactForceDerivative.zy()*w0;
                                lcCoeffs[bI][segI][nbI].first()(2, 2) -=
                                    curContactForceDerivative.zz()*w0;

                                //- upper
                                lcCoeffs[bI][segI][nbI].second()(0, 0) -=
                                    curContactForceDerivative.xx()*w1;
                                lcCoeffs[bI][segI][nbI].second()(0, 1) -=
                                    curContactForceDerivative.xy()*w1;
                                lcCoeffs[bI][segI][nbI].second()(0, 2) -=
                                    curContactForceDerivative.xz()*w1;

                                lcCoeffs[bI][segI][nbI].second()(1, 0) -=
                                    curContactForceDerivative.yx()*w1;
                                lcCoeffs[bI][segI][nbI].second()(1, 1) -=
                                    curContactForceDerivative.yy()*w1;
                                lcCoeffs[bI][segI][nbI].second()(1, 2) -=
                                    curContactForceDerivative.yz()*w1;

                                lcCoeffs[bI][segI][nbI].second()(2, 0) -=
                                    curContactForceDerivative.zx()*w1;
                                lcCoeffs[bI][segI][nbI].second()(2, 1) -=
                                    curContactForceDerivative.zy()*w1;
                                lcCoeffs[bI][segI][nbI].second()(2, 2) -=
                                    curContactForceDerivative.zz()*w1;


                                //   // Contact force derivative (Theta)

                                //   //- lower
                                //   lcCoeffs[bI][segI][nbI].first()(0, 3) +=
                                    //   curCircumContactForceThetaDerivative.xx()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(0, 4) +=
                                    //   curCircumContactForceThetaDerivative.xy()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(0, 5) +=
                                    //   curCircumContactForceThetaDerivative.xz()
                                   //   *L()[cellID]*w0;

                                //   lcCoeffs[bI][segI][nbI].first()(1, 3) +=
                                    //   curCircumContactForceThetaDerivative.yx()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(1, 4) +=
                                    //   curCircumContactForceThetaDerivative.yy()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(1, 5) +=
                                    //   curCircumContactForceThetaDerivative.yz()
                                   //   *L()[cellID]*w0;

                                //   lcCoeffs[bI][segI][nbI].first()(2, 3) +=
                                    //   curCircumContactForceThetaDerivative.zx()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(2, 4) +=
                                    //   curCircumContactForceThetaDerivative.zy()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(2, 5) +=
                                    //   curCircumContactForceThetaDerivative.zz()
                                   //   *L()[cellID]*w0;

                                //   //- upper
                                //   lcCoeffs[bI][segI][nbI].second()(0, 3) +=
                                    //   curCircumContactForceThetaDerivative.xx()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(0, 4) +=
                                    //   curCircumContactForceThetaDerivative.xy()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(0, 5) +=
                                    //   curCircumContactForceThetaDerivative.xz()
                                   //   *L()[cellID]*w1;

                                //   lcCoeffs[bI][segI][nbI].second()(1, 3) +=
                                    //   curCircumContactForceThetaDerivative.yx()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(1, 4) +=
                                    //   curCircumContactForceThetaDerivative.yy()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(1, 5) +=
                                    //   curCircumContactForceThetaDerivative.yz()
                                   //   *L()[cellID]*w1;

                                //   lcCoeffs[bI][segI][nbI].second()(2, 3) +=
                                    //   curCircumContactForceThetaDerivative.zx()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(2, 4) +=
                                    //   curCircumContactForceThetaDerivative.zy()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(2, 5) +=
                                    //   curCircumContactForceThetaDerivative.zz()
                                   //   *L()[cellID]*w1;

                                // Contact moment derivative

                                //   //- lower
                                //   lcCoeffs[bI][segI][nbI].first()(3, 3) +=
                                    //   curContactMomentDerivative.xx()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(3, 4) +=
                                    //   curContactMomentDerivative.xy()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(3, 5) +=
                                    //   curContactMomentDerivative.xz()
                                   //   *L()[cellID]*w0;

                                //   lcCoeffs[bI][segI][nbI].first()(4, 3) +=
                                    //   curContactMomentDerivative.yx()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(4, 4) +=
                                    //   curContactMomentDerivative.yy()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(4, 5) +=
                                    //   curContactMomentDerivative.yz()
                                   //   *L()[cellID]*w0;

                                //   lcCoeffs[bI][segI][nbI].first()(5, 3) +=
                                    //   curContactMomentDerivative.zx()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(5, 4) +=
                                    //   curContactMomentDerivative.zy()
                                   //   *L()[cellID]*w0;
                                //   lcCoeffs[bI][segI][nbI].first()(5, 5) +=
                                    //   curContactMomentDerivative.zz()
                                   //   *L()[cellID]*w0;

                                //   //- upper
                                //   lcCoeffs[bI][segI][nbI].second()(3, 3) +=
                                    //   curContactMomentDerivative.xx()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(3, 4) +=
                                    //   curContactMomentDerivative.xy()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(3, 5) +=
                                    //   curContactMomentDerivative.xz()
                                   //   *L()[cellID]*w1;

                                //   lcCoeffs[bI][segI][nbI].second()(4, 3) +=
                                    //   curContactMomentDerivative.yx()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(4, 4) +=
                                    //   curContactMomentDerivative.yy()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(4, 5) +=
                                    //   curContactMomentDerivative.yz()
                                   //   *L()[cellID]*w1;

                                //   lcCoeffs[bI][segI][nbI].second()(5, 3) +=
                                    //   curContactMomentDerivative.zx()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(5, 4) +=
                                    //   curContactMomentDerivative.zy()
                                   //   *L()[cellID]*w1;
                                //   lcCoeffs[bI][segI][nbI].second()(5, 5) +=
                                    //   curContactMomentDerivative.zz()
                                   //   *L()[cellID]*w1;
                            }
                        }
                    }
                }
	    }
	    start += nSeg;
	}
    }
}



void coupledTotalLagNewtonRaphsonBeam::applyPointContact
(
    Field<scalarSquareMatrix>& diag,
    Field<scalarSquareMatrix>& lower,
    Field<scalarSquareMatrix>& upper,
    Field<scalarRectangularMatrix>& source
    // multibeamFvBlockMatrix& eqn
)
{
    if
    (
        (mesh().cellZones().size() > 1)
     && !Pstream::parRun()
    )
    {
    // Info << "Inside applyPointContact() in applyContacts.C \n" << endl;
    // // Get matrix diagonal
    // tensor6Field& diag = eqn.diag().asSquare();

    // // Get matrix upper coeffs
    // tensor6Field& upper = eqn.upper().asSquare();

    // // Get matrix lower coeffs
    // tensor6Field& lower = eqn.lower().asSquare();

    // // Grab source
    // vector6Field& source = eqn.source();

    // Grab beam coupling coeffs (foam-extend version)
    //    tensor6PairListList& pcCoeffs = eqn.pointContactCoeffs();

    // SB (2026 ESI) - A storage container to collect lcCoeffs for ESI version
    // IMPORTANT - However, these coefficients are simply collected but not actually
    // added to the block matrix of beam solver -
    // The lcCoeffs (and pcCoeffs) were passed to multiBeamFvMatrix class where
    // contributions were added. Since we do not use the class, the contribution of
    // the coeffs is simply being lost at this point.
    pcCoeffPairListList pcCoeffs(contact().pointContacts().size());

    const pcCoeffPair pcZero
    (
        scalarRectangularMatrix(6, 3, Zero),
        scalarRectangularMatrix(6, 3, Zero)
    );

    forAll(pcCoeffs, pcI)
    {
        pcCoeffs[pcI].setSize(2, pcZero);
    }

    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    const surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

    const surfaceScalarField& dc = mesh().deltaCoeffs();

    // label start = 0;
    forAll(contact().pointContacts(), pcI)
    {
	//   Info << "\nInside contact().pointContacts() in applyContacts.C" << endl;
        labelList bI(2, -1);
        labelList globalSegI(2, -1);
        scalarField zeta(2, -2);
        vectorList contactForce(2, vector::zero);
        vectorList tangContactForce(2, vector::zero);

	//   volVectorField Wprev = W_.prevIter();

	//   Info << "W:\n " << Wprev << "\n\nW: " << W_ << endl;


        //Owner
        {
            bI[0] = contact().pointContacts()[pcI].firstBeam();
            label segI = contact().pointContacts()[pcI].firstBeamSegment();
            zeta[0] = contact().pointContacts()[pcI].firstBeamZeta();

            // label start = 0;
            // label i = 0;
            // while(i < bI[0])
            // {
            //     start += contact().splines()[i].nSegments();
            //     i++;
            // }
            // globalSegI[0] = start + segI;

            globalSegI[0] = whichCell(bI[0], segI);

	    // SB (2023 May): changed the explicit contact contribution
            contactForce[0] =
                contact().pointContacts()[pcI].normalContactForce();
	    tangContactForce[0] =
                contact().pointContacts()[pcI].firstBeamTangContactForce()
              - contact().pointContacts()[pcI].secondBeamTangContactForce();

	   //   Info << "segI " << segI << endl;

        }

        // Neighbour
        {
            bI[1] = contact().pointContacts()[pcI].secondBeam();
            label neiSegI = contact().pointContacts()[pcI].secondBeamSegment();
            zeta[1] = contact().pointContacts()[pcI].secondBeamZeta();

            // label neiStart = 0;
            // label i = 0;
            // while(i < bI[1])
            // {
            //     neiStart += contact().splines()[i].nSegments();
            //     i++;
            // }
            // globalSegI[1] = neiStart + neiSegI;

            globalSegI[1] = whichCell(bI[1], neiSegI);


	    // SB (2023 May): changed the explicit contact contribution
            contactForce[1] = -contactForce[0];
		//   - contact().pointContacts()[pcI].normalContactForce()
	    tangContactForce[1] =
		  contact().pointContacts()[pcI].secondBeamTangContactForce()
		- contact().pointContacts()[pcI].firstBeamTangContactForce();


        }


        forAll(globalSegI, sI)
        {
            vector DR = vector::zero;

            if (zeta[sI] > 0)
            {
                // label faceID = findIndex(own, globalSegI[sI]);
                label faceID = own.find(globalSegI[sI]);

                //		Info << "faceID own: " << faceID << endl;
                if (faceID == -1) // last cell
                {
                    // const unallocLabelList& faceCells =
                    const labelUList& faceCells =                        
                        mesh().boundary()[endPatchIndex(bI[sI])].faceCells();

                    // label bFaceID = findIndex(faceCells, globalSegI[sI]);
                    label bFaceID = faceCells.find(globalSegI[sI]);

                    DR = zeta[sI]
                       *dRdS.boundaryField()[endPatchIndex(bI[sI])][bFaceID]
                       /dc.boundaryField()[endPatchIndex(bI[sI])][bFaceID];
                }
                else
                {
                    DR = 0.5*zeta[sI]*dRdS.internalField()[faceID]
                       /dc.internalField()[faceID];

                }
            }
            else
            {
                // label faceID = findIndex(nei, globalSegI[sI]);
                label faceID = nei.find(globalSegI[sI]);                
                //		Info << "faceID nei: " << faceID << endl;
                if (faceID == -1) // first cell
                {
                    // const unallocLabelList& faceCells =
                    const labelUList& faceCells =                    
                        mesh().boundary()[startPatchIndex(bI[sI])].faceCells();

                    // label bFaceID = findIndex(faceCells, globalSegI[sI]);
                    label bFaceID = faceCells.find(globalSegI[sI]);                    

                    DR = zeta[sI]
                       *dRdS.boundaryField()[startPatchIndex(bI[sI])][bFaceID]
                       /dc.boundaryField()[startPatchIndex(bI[sI])][bFaceID];
                }
                else
                {
                    DR = 0.5*zeta[sI]*dRdS.internalField()[faceID]
                       /dc.internalField()[faceID];

                }
            }

            vector F0 = contactForce[sI] + tangContactForce[sI];

            vector M0 =  (spinTensor(DR) & F0);




            // label index = 6*globalSegI[sI];
            // b(index++) -= F0.x();
            // b(index++) -= F0.y();
            // b(index++) -= F0.z();
            // b(index++) -= M0.x();
            // b(index++) -= M0.y();
            // b(index++) -= M0.z();

            // W equation
            source[globalSegI[sI]](0, 0) -= F0.x();
            source[globalSegI[sI]](1, 0) -= F0.y();
            source[globalSegI[sI]](2, 0) -= F0.z();


            // Theta equation
            source[globalSegI[sI]](3, 0) -= M0.x();
            source[globalSegI[sI]](4, 0) -= M0.y();
            source[globalSegI[sI]](5, 0) -= M0.z();
        }
    }

    //Info << "Explicit part of point force calculation done " << endl;

    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A;

    const cellList& cells = mesh().cells();
    // const labelList& own = mesh().owner();

    // Implicit part of the point contact force
    forAll(contact().pointContacts(), pcI)
    {
        labelList bI(2, -1);
        labelList nbI(2, -1);
        labelList segI(2, -1);
        labelList neiSegI(2, -1);
        labelList globalSegI(2, -1);
        labelList globalLowerSegI(2, -1);
        labelList globalUpperSegI(2, -1);
        scalarField weight(2, 0);
        labelList globalNeiSegI(2, -1);
        labelList globalLowerNeiSegI(2, -1);
        labelList globalUpperNeiSegI(2, -1);
        scalarField neiWeight(2, 0);
        scalarField zeta(2, -2);
        scalarField neiZeta(2, -2);
        scalarField contactDistance(2, 0);
        //   vectorList contactForce(2, vector::zero);
        tensorList contactForceDerivative(2, tensor::zero);
	//SB (2023 May): Added contact force dn component
	tensorList contactForceDn(2, tensor::zero);

        vectorList tangentialContactForce(2, vector::zero);
        tensorList tangContactForceDerivative(2, tensor::zero);
	tensorList tangContactForceDt(2,tensor::zero);

	tensorList couplingTangentContactForceCoeffs(2,tensor::zero);


        //Owner
        {
            bI[0] = contact().pointContacts()[pcI].firstBeam();
            segI[0] = contact().pointContacts()[pcI].firstBeamSegment();
            zeta[0] = contact().pointContacts()[pcI].firstBeamZeta();
            globalSegI[0] = whichCell(bI[0], segI[0]);
            globalLowerSegI[0] =
                whichCell
                (
                    bI[0],
                    contact().pointContacts()[pcI].firstBeamLowerSegment()
                );
            globalUpperSegI[0] =
                whichCell
                (
                    bI[0],
                    contact().pointContacts()[pcI].firstBeamUpperSegment()
                );

            weight[0] = contact().pointContacts()[pcI].firstBeamWeight();

            contactDistance[0] =
                contact().pointContacts()[pcI].delta();
            //   contactForce[0] =
                //   contact().pointContacts()[pcI].normalContactForce();
            contactForceDerivative[0] =
                contact().pointContacts()[pcI].normalContactForceDerivative();

	    //SB (2023 May): Added contact force dn component
	    contactForceDn[0] = contact().pointContacts()[pcI].normalContactForceDn();

            tangentialContactForce[0] =
                contact().pointContacts()[pcI].firstBeamTangContactForce();

            //   tangContactForceDerivative[0] =
                //   contact().pointContacts()[pcI]
               //   .firstBeamTangContactForceDerivative();
              //   - contact().pointContacts()[pcI]
               //   .secondBeamTangContactForceDerivative();

	    /*--SB (2023 May): changed the frictional contact contributions---*/

	    tangContactForceDerivative[0] =
		contact().pointContacts()[pcI].firstBeamTangContactForceDerivative()
	      + contact().pointContacts()[pcI].couplingFirstBeamTangContactForceDerivative();


	    tangContactForceDt[0] =
		contact().pointContacts()[pcI].firstBeamTangContactForceDt()
	      + contact().pointContacts()[pcI].couplingFirstBeamTangContactForceDt();

	    couplingTangentContactForceCoeffs[0] =
		contact().pointContacts()[pcI].secondBeamTangContactForceDerivative()
	      + contact().pointContacts()[pcI].secondBeamTangContactForceDt()
	      + contact().pointContacts()[pcI].couplingSecondBeamTangContactForceDerivative()
	      + contact().pointContacts()[pcI].couplingSecondBeamTangContactForceDt();

	    /*---*/
        }

        // Neighbour
        {
            bI[1] = contact().pointContacts()[pcI].secondBeam();
            segI[1] = contact().pointContacts()[pcI].secondBeamSegment();
            zeta[1] = contact().pointContacts()[pcI].secondBeamZeta();
            globalSegI[1] = whichCell(bI[1], segI[1]);
            globalLowerSegI[1] =
                whichCell
                (
                    bI[1],
                    contact().pointContacts()[pcI].secondBeamLowerSegment()
                );
            globalUpperSegI[1] =
                whichCell
                (
                    bI[1],
                    contact().pointContacts()[pcI].secondBeamUpperSegment()
                );

            weight[1] = contact().pointContacts()[pcI].secondBeamWeight();

            contactDistance[1] = contactDistance[0];
            //   contactForce[1] = -contactForce[0];

            contactForceDerivative[1] = contactForceDerivative[0];

	    //SB (2023 May): Added contact force dn component
	    contactForceDn[1] = contactForceDn[0];

            tangentialContactForce[1] =
                contact().pointContacts()[pcI].secondBeamTangContactForce();

            //   tangContactForceDerivative[1] =
                //   tangContactForceDerivative[0];

	    /*--SB (2023 May): changed the frictional contact contributions---*/

	    tangContactForceDerivative[1] =
		contact().pointContacts()[pcI].secondBeamTangContactForceDerivative()
	      + contact().pointContacts()[pcI].couplingSecondBeamTangContactForceDerivative();


	    tangContactForceDt[1] =
		contact().pointContacts()[pcI].secondBeamTangContactForceDt()
	      + contact().pointContacts()[pcI].couplingSecondBeamTangContactForceDt();

	    //   tangContactForceDerivative[1] =
		//   contact().pointContacts()[pcI].secondBeamTangContactForceDerivative();


	    //   tangContactForceDt[1] = //tangContactForceDt[0];
		//   contact().pointContacts()[pcI].secondBeamTangContactForceDt();

	    couplingTangentContactForceCoeffs[1] =
		contact().pointContacts()[pcI].firstBeamTangContactForceDerivative()
	      + contact().pointContacts()[pcI].firstBeamTangContactForceDt()
	      + contact().pointContacts()[pcI].couplingFirstBeamTangContactForceDerivative()
	      + contact().pointContacts()[pcI].couplingFirstBeamTangContactForceDt();

	    /*---*/
        }

        nbI[0] = bI[1];
        nbI[1] = bI[0];

        neiSegI[0] = segI[1];
        neiSegI[1] = segI[0];

        globalNeiSegI[0] = globalSegI[1];
        globalNeiSegI[1] = globalSegI[0];

        globalLowerNeiSegI[0] = globalLowerSegI[1];
        globalLowerNeiSegI[1] = globalLowerSegI[0];

        globalUpperNeiSegI[0] = globalUpperSegI[1];
        globalUpperNeiSegI[1] = globalUpperSegI[0];

        neiWeight[0] = weight[1];
        neiWeight[1] = weight[0];

        neiZeta[0] = zeta[1];
        neiZeta[1] = zeta[0];

        forAll(globalSegI, sI)
        {
	    //   Info << "sI: " << sI << endl;
            vector DR = vector::zero;
            if (zeta[sI] > 0)
            {
                // label faceID = findIndex(own, globalSegI[sI]);
                label faceID = own.find(globalSegI[sI]);
                if (faceID == -1) // last cell
                {
                    // const unallocLabelList& faceCells =
                    const labelUList& faceCells =                        
                        mesh().boundary()[endPatchIndex(bI[sI])].faceCells();

                    // label bFaceID = findIndex(faceCells, globalSegI[sI]);
                    label bFaceID = faceCells.find(globalSegI[sI]);
                    
                    DR = zeta[sI]
                       *dRdS.boundaryField()[endPatchIndex(bI[sI])][bFaceID]
                       /dc.boundaryField()[endPatchIndex(bI[sI])][bFaceID];
                }
                else
                {
                    DR = 0.5*zeta[sI]*dRdS.internalField()[faceID]
                       /dc.internalField()[faceID];
                }
            }
            else
            {
                // label faceID = findIndex(nei, globalSegI[sI]);
                label faceID = nei.find(globalSegI[sI]);
                if (faceID == -1) // first cell
                {
                    // const unallocLabelList& faceCells =
                    const labelUList& faceCells =                        
                        mesh().boundary()[startPatchIndex(bI[sI])].faceCells();

                    // label bFaceID =
                    //     findIndex(faceCells, globalSegI[sI]);
                    label bFaceID = faceCells.find(globalSegI[sI]);

                    DR = zeta[sI]
                       *dRdS.boundaryField()[startPatchIndex(bI[sI])][bFaceID]
                       /dc.boundaryField()[startPatchIndex(bI[sI])][bFaceID];
                }
                else
                {
                    DR = 0.5*zeta[sI]*dRdS.internalField()[faceID]
                       /dc.internalField()[faceID];
                }
            }

            tensor DMCoeff =
	    (
		spinTensor(DR)
		& (contactForceDerivative[sI] + contactForceDn[sI])
	    );

	    //   tensor DMCoeff = tensor::zero;


            scalar wDiag = weight[sI];
            label diagCell = globalLowerSegI[sI];


            if (globalSegI[sI] != globalLowerSegI[sI])
            {
                wDiag = 1.0 - wDiag;
                diagCell = globalUpperSegI[sI];

            }
            label sharedFaceIndex =
                sharedFace
                (
                    cells[globalLowerSegI[sI]],
                    cells[globalUpperSegI[sI]]
                );

            //-- wDiag
            //- Force
            diag[diagCell](0, 0) +=
                wDiag*contactForceDerivative[sI].xx()
              + wDiag*tangContactForceDerivative[sI].xx()
	      + wDiag*contactForceDn[sI].xx()
	      + wDiag*tangContactForceDt[sI].xx();

            diag[diagCell](0, 1) +=
                wDiag*contactForceDerivative[sI].xy()
              + wDiag*tangContactForceDerivative[sI].xy()
	      + wDiag*contactForceDn[sI].xy()
              + wDiag*tangContactForceDt[sI].xy();

            diag[diagCell](0, 2) +=
                wDiag*contactForceDerivative[sI].xz()
              + wDiag*tangContactForceDerivative[sI].xz()
	      + wDiag*contactForceDn[sI].xz()
              + wDiag*tangContactForceDt[sI].xz();

            diag[diagCell](1, 0) +=
                wDiag*contactForceDerivative[sI].yx()
              + wDiag*tangContactForceDerivative[sI].yx()
	      + wDiag*contactForceDn[sI].yx()
	      + wDiag*tangContactForceDt[sI].yx();

            diag[diagCell](1, 1) +=
                wDiag*contactForceDerivative[sI].yy()
              + wDiag*tangContactForceDerivative[sI].yy()
	      + wDiag*contactForceDn[sI].yy()
	      + wDiag*tangContactForceDt[sI].yy();

            diag[diagCell](1, 2) +=
                wDiag*contactForceDerivative[sI].yz()
              + wDiag*tangContactForceDerivative[sI].yz()
	      + wDiag*contactForceDn[sI].yz()
	      + wDiag*tangContactForceDt[sI].yz();


            diag[diagCell](2, 0) +=
                wDiag*contactForceDerivative[sI].zx()
              + wDiag*tangContactForceDerivative[sI].zx()
	      + wDiag*contactForceDn[sI].zx()
	      + wDiag*tangContactForceDt[sI].zx();

            diag[diagCell](2, 1) +=
                wDiag*contactForceDerivative[sI].zy()
              + wDiag*tangContactForceDerivative[sI].zy()
              + wDiag*contactForceDn[sI].zy()
	      + wDiag*tangContactForceDt[sI].zy();

            diag[diagCell](2, 2) +=
                wDiag*contactForceDerivative[sI].zz()
              + wDiag*tangContactForceDerivative[sI].zz()
              + wDiag*contactForceDn[sI].zz()
	      + wDiag*tangContactForceDt[sI].zz();


            //- Moment
            diag[diagCell](3, 0) += wDiag*DMCoeff.xx();
            diag[diagCell](3, 1) += wDiag*DMCoeff.xy();
            diag[diagCell](3, 2) += wDiag*DMCoeff.xz();

            diag[diagCell](4, 0) += wDiag*DMCoeff.yx();
            diag[diagCell](4, 1) += wDiag*DMCoeff.yy();
            diag[diagCell](4, 2) += wDiag*DMCoeff.yz();

            diag[diagCell](5, 0) += wDiag*DMCoeff.zx();
            diag[diagCell](5, 1) += wDiag*DMCoeff.zy();
            diag[diagCell](5, 2) += wDiag*DMCoeff.zz();

            //-- w1
            if (sharedFaceIndex != -1)
            {
                scalarSquareMatrix* curCoeffPtr = &(upper[sharedFaceIndex]);
                if (own[sharedFaceIndex] != globalSegI[sI])
                {
                    curCoeffPtr = &(lower[sharedFaceIndex]);
                }
                scalarSquareMatrix& curCoeff = *curCoeffPtr;

                //- Force
                curCoeff(0, 0) +=
                    (1.0-wDiag)*contactForceDerivative[sI].xx()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].xx()
                  + (1.0-wDiag)*contactForceDn[sI].xx()
		  + (1.0-wDiag)*tangContactForceDt[sI].xx();

                curCoeff(0, 1) +=
                    (1.0-wDiag)*contactForceDerivative[sI].xy()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].xy()
                  + (1.0-wDiag)*contactForceDn[sI].xy()
		  + (1.0-wDiag)*tangContactForceDt[sI].xy();

                curCoeff(0, 2) +=
                    (1.0-wDiag)*contactForceDerivative[sI].xz()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].xz()
                  + (1.0-wDiag)*contactForceDn[sI].xz()
		  + (1.0-wDiag)*tangContactForceDt[sI].xz();


                curCoeff(1, 0) +=
                    (1.0-wDiag)*contactForceDerivative[sI].yx()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].yx()
                  + (1.0-wDiag)*contactForceDn[sI].yx()
		  + (1.0-wDiag)*tangContactForceDt[sI].yx();

                curCoeff(1, 1) +=
                    (1.0-wDiag)*contactForceDerivative[sI].yy()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].yy()
                  + (1.0-wDiag)*contactForceDn[sI].yy()
		  + (1.0-wDiag)*tangContactForceDt[sI].yy();

                curCoeff(1, 2) +=
                    (1.0-wDiag)*contactForceDerivative[sI].yz()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].yz()
                  + (1.0-wDiag)*contactForceDn[sI].yz()
		  + (1.0-wDiag)*tangContactForceDt[sI].yz();


                curCoeff(2, 0) +=
                    (1.0-wDiag)*contactForceDerivative[sI].zx()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].zx()
                  + (1.0-wDiag)*contactForceDn[sI].zx()
		  + (1.0-wDiag)*tangContactForceDt[sI].zx();

                curCoeff(2, 1) +=
                    (1.0-wDiag)*contactForceDerivative[sI].zy()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].zy()
                  + (1.0-wDiag)*contactForceDn[sI].zy()
		  + (1.0-wDiag)*tangContactForceDt[sI].zy();

                curCoeff(2, 2) +=
                    (1.0-wDiag)*contactForceDerivative[sI].zz()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].zz()
                  + (1.0-wDiag)*contactForceDn[sI].zz()
		  + (1.0-wDiag)*tangContactForceDt[sI].zz();


                //- Moment
                curCoeff(3, 0) +=
                    (1.0-wDiag)*DMCoeff.xx();
                curCoeff(3, 1) +=
                    (1.0-wDiag)*DMCoeff.xy();
                curCoeff(3, 2) +=
                    (1.0-wDiag)*DMCoeff.xz();

                curCoeff(4, 0) +=
                    (1.0-wDiag)*DMCoeff.yx();
                curCoeff(4, 1) +=
                    (1.0-wDiag)*DMCoeff.yy();
                curCoeff(4, 2) +=
                    (1.0-wDiag)*DMCoeff.yz();

                curCoeff(5, 0) +=
                    (1.0-wDiag)*DMCoeff.zx();
                curCoeff(5, 1) +=
                    (1.0-wDiag)*DMCoeff.zy();
                curCoeff(5, 2) +=
                    (1.0-wDiag)*DMCoeff.zz();
            }



            // Off-diagonal
            // label globalNeiSeg0 = globalLowerNeiSegI[sI];
            // label globalNeiSeg1 = globalUpperNeiSegI[sI];
            // label globalNeiSeg0 = globalNeiSegI[sI];
            // label globalNeiSeg1 = globalNeiSegI[sI];

            scalar w0 = neiWeight[sI];
            scalar w1 = 1.0 - w0;


            // label nNeiSegments =
            //     contact().splines()[nbI[sI]].nSegments();
            // scalarField neiL =
            //     contact().splines()[nbI[sI]].segLengths();

            // if
            // (
            //     (neiSegI[sI] > 0)
            //  && (neiSegI[sI] < (nNeiSegments-1))
            // )
            // {
            //     if (neiZeta[sI] >= 0)
            //     {
            //         globalNeiSeg1 = globalNeiSegI[sI] + 1;
            //         scalar l = 0.5*neiL[neiSegI[sI]]
            //           + 0.5*neiL[neiSegI[sI]+1];
            //         scalar l0 = neiZeta[sI]*neiL[neiSegI[sI]]/2;
            //         w0 = (l-l0)/l;
            //         w1 = 1-w0;
            //     }
            //     else
            //     {
            //         globalNeiSeg0 = globalNeiSegI[sI] - 1;
            //         scalar l = 0.5*neiL[neiSegI[sI]-1]
            //           + 0.5*neiL[neiSegI[sI]];
            //         scalar l1 = -neiZeta[sI]*neiL[neiSegI[sI]]/2;
            //         w1 = (l-l1)/l;
            //         w0 = 1-w1;
            //     }
            // }

            // label IN0 = 6*globalNeiSeg0;
            // label IN1 = 6*globalNeiSeg1;

            //-- IN0
            //- Force
            pcCoeffs[pcI][sI].first()(0,0) -=
                w0*contactForceDerivative[sI].xx()
              + w0*contactForceDn[sI].xx()
	      + w0*couplingTangentContactForceCoeffs[sI].xx();
            pcCoeffs[pcI][sI].first()(0,1) -=
                w0*contactForceDerivative[sI].xy()
              + w0*contactForceDn[sI].xy()
	      + w0*couplingTangentContactForceCoeffs[sI].xy();
            pcCoeffs[pcI][sI].first()(0,2) -=
                w0*contactForceDerivative[sI].xz()
              + w0*contactForceDn[sI].xz()
	      + w0*couplingTangentContactForceCoeffs[sI].xz();

            pcCoeffs[pcI][sI].first()(1,0) -=
                w0*contactForceDerivative[sI].yx()
              + w0*contactForceDn[sI].yx()
	      + w0*couplingTangentContactForceCoeffs[sI].yx();
            pcCoeffs[pcI][sI].first()(1,1) -=
                w0*contactForceDerivative[sI].yy()
              + w0*contactForceDn[sI].yy()
	      + w0*couplingTangentContactForceCoeffs[sI].yy();
            pcCoeffs[pcI][sI].first()(1,2) -=
                w0*contactForceDerivative[sI].yz()
              + w0*contactForceDn[sI].yz()
	      + w0*couplingTangentContactForceCoeffs[sI].yz();

            pcCoeffs[pcI][sI].first()(2,0) -=
                w0*contactForceDerivative[sI].zx()
              + w0*contactForceDn[sI].zx()
	      + w0*couplingTangentContactForceCoeffs[sI].zx();
            pcCoeffs[pcI][sI].first()(2,1) -=
                w0*contactForceDerivative[sI].zy()
              + w0*contactForceDn[sI].zy()
	      + w0*couplingTangentContactForceCoeffs[sI].zy();
            pcCoeffs[pcI][sI].first()(2,2) -=
                w0*contactForceDerivative[sI].zz()
              + w0*contactForceDn[sI].zz()
	      + w0*couplingTangentContactForceCoeffs[sI].zz();

	    //Info << "pcCoeffs force: " << pcCoeffs << endl;
            //- Moment
            pcCoeffs[pcI][sI].first()(3,0) -=
                w0*DMCoeff.xx();
            pcCoeffs[pcI][sI].first()(3,1) -=
                w0*DMCoeff.xy();
            pcCoeffs[pcI][sI].first()(3,2) -=
                w0*DMCoeff.xz();

            pcCoeffs[pcI][sI].first()(4,0) -=
                w0*DMCoeff.yx();
            pcCoeffs[pcI][sI].first()(4,1) -=
                w0*DMCoeff.yy();
            pcCoeffs[pcI][sI].first()(4,2) -=
                w0*DMCoeff.yz();

            pcCoeffs[pcI][sI].first()(5,0) -=
                w0*DMCoeff.zx();
            pcCoeffs[pcI][sI].first()(5,1) -=
                w0*DMCoeff.zy();
            pcCoeffs[pcI][sI].first()(5,2) -=
                w0*DMCoeff.zz();
	    //Info << "pcCoeffs moment: " << pcCoeffs << endl;
            //-- IN1
            //- Force
            pcCoeffs[pcI][sI].second()(0,0) -=
                w1*contactForceDerivative[sI].xx()
              + w1*contactForceDn[sI].xx()
	      + w1*couplingTangentContactForceCoeffs[sI].xx();
            pcCoeffs[pcI][sI].second()(0,1) -=
                w1*contactForceDerivative[sI].xy()
              + w1*contactForceDn[sI].xy()
	      + w1*couplingTangentContactForceCoeffs[sI].xy();
            pcCoeffs[pcI][sI].second()(0,2) -=
                w1*contactForceDerivative[sI].xz()
              + w1*contactForceDn[sI].xz()
	      + w1*couplingTangentContactForceCoeffs[sI].xz();

            pcCoeffs[pcI][sI].second()(1,0) -=
                w1*contactForceDerivative[sI].yx()
              + w1*contactForceDn[sI].yx()
	      + w1*couplingTangentContactForceCoeffs[sI].yx();
            pcCoeffs[pcI][sI].second()(1,1) -=
                w1*contactForceDerivative[sI].yy()
              + w1*contactForceDn[sI].yy()
	      + w1*couplingTangentContactForceCoeffs[sI].yy();
            pcCoeffs[pcI][sI].second()(1,2) -=
                w1*contactForceDerivative[sI].yz()
              + w1*contactForceDn[sI].yz()
	      + w1*couplingTangentContactForceCoeffs[sI].yz();

            pcCoeffs[pcI][sI].second()(2,0) -=
                w1*contactForceDerivative[sI].zx()
              + w1*contactForceDn[sI].zx()
	      + w1*couplingTangentContactForceCoeffs[sI].zx();
            pcCoeffs[pcI][sI].second()(2,1) -=
                w1*contactForceDerivative[sI].zy()
              + w1*contactForceDn[sI].zy()
	      + w1*couplingTangentContactForceCoeffs[sI].zy();
            pcCoeffs[pcI][sI].second()(2,2) -=
                w1*contactForceDerivative[sI].zz()
              + w1*contactForceDn[sI].zz()
	      + w1*couplingTangentContactForceCoeffs[sI].zz();

            //- Moment
            pcCoeffs[pcI][sI].second()(3,0) -=
                w1*DMCoeff.xx();
            pcCoeffs[pcI][sI].second()(3,1) -=
                w1*DMCoeff.xy();
            pcCoeffs[pcI][sI].second()(3,2) -=
                w1*DMCoeff.xz();

            pcCoeffs[pcI][sI].second()(4,0) -=
                w1*DMCoeff.yx();
            pcCoeffs[pcI][sI].second()(4,1) -=
                w1*DMCoeff.yy();
            pcCoeffs[pcI][sI].second()(4,2) -=
                w1*DMCoeff.yz();

            pcCoeffs[pcI][sI].second()(5,0) -=
                w1*DMCoeff.zx();
            pcCoeffs[pcI][sI].second()(5,1) -=
                w1*DMCoeff.zy();
            pcCoeffs[pcI][sI].second()(5,2) -=
                w1*DMCoeff.zz();

	    //Info << "pcCoeffs: " << pcCoeffs << endl;
            // // Tangential force direction corrector
            // if (false)
            // {
            //     scalar delta;
            //     label ID = 6*globalSegI[sI];
            //     label ID1 = 6*globalSeg1;
            //     label IND = globalNeiSegI[sI];

            //     if (globalSeg1 > globalSegI[sI])
            //     {
            //         delta = (0.5*L()[segI[sI]] + 0.5*L()[segI[sI]+1]);
            //     }
            //     else
            //     {
            //         delta = -(0.5*L()[segI[sI]-1] + 0.5*L()[segI[sI]]);
            //     }

            //     //
            //     A.coeffRef(ID,ID) -=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(ID+1,ID+1) -=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(ID+2,ID+2) -=
            //         mag(tangentialContactForce[sI])/delta;

            //     A.coeffRef(ID,ID1) +=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(ID+1,ID1+1) +=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(ID+2,ID1+2) +=
            //         mag(tangentialContactForce[sI])/delta;

            //     //
            //     A.coeffRef(IND,ID) +=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(IND+1,ID+1) +=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(IND+2,ID+2) +=
            //         mag(tangentialContactForce[sI])/delta;

            //     A.coeffRef(IND,ID1) -=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(IND+1,ID1+1) -=
            //         mag(tangentialContactForce[sI])/delta;
            //     A.coeffRef(IND+2,ID1+2) -=
            //         mag(tangentialContactForce[sI])/delta;
            // }
        }
    }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
