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
#include "lagrangeMultipliers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace beamModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void coupledTotalLagNewtonRaphsonBeam::applyLineContact
(
    multibeamFvBlockMatrix& eqn
)
{
    if (mesh().cellZones().size() > 1)
    {
    	// Info << "Inside applyLineContact() in applyContacts.C \n" << endl;
        label nBeams = contact().splines().size();
        
        // Get matrix diagonal
        tensor6Field& diag = eqn.diag().asSquare();
        
        // Grab source
        vector6Field& source = eqn.source();

        // Grab beam coupling coeffs
        PtrList<tensor6PairListList>& lcCoeffs = eqn.lineContactCoeffs();

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

            for
            (
                label segI=0;
                segI<contact().splines()[bI].nSegments();
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
                        }
                    }
                }
            }

            forAll(curq, segI)
            {
                label globalSegIndex = start + segI;

                label cellID =
                    localCellIndex(globalSegIndex);

                if (cellID != -1)
                {
                    // W equation
                    source[cellID](0) -= curq[segI].x()*L()[cellID];
                    source[cellID](1) -= curq[segI].y()*L()[cellID];
                    source[cellID](2) -= curq[segI].z()*L()[cellID];

                    // Theta equation
                    source[cellID](3) -= curm[segI].x()*L()[cellID];
                    source[cellID](4) -= curm[segI].y()*L()[cellID];
                    source[cellID](5) -= curm[segI].z()*L()[cellID];
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
                            const label neiSegI =
                                curLineContacts[segI][nbI].secondBeamSegment();

                            const vector& curContactForce =
                                curLineContacts[segI][nbI].normalContactForce();

                            if (mag(curContactForce) > SMALL)
                            {
                                const scalar neiSegParam =
                                    curLineContacts[segI][nbI].secondBeamZeta();

                                // const label neiLowerSegI =
                                //     curLineContacts[segI][nbI]
                                //    .secondBeamLowerSegment();
                                // const label neiUpperSegI =
                                //     curLineContacts[segI][nbI]
                                //    .secondBeamUpperSegment();

                                const scalar& curContactDistance =
                                    curLineContacts[segI][nbI].delta();

                                tensor curContactForceDerivative =
                                    curLineContacts[segI][nbI]
                                   .normalContactForceDerivative()
                                   *L()[cellID];

                                const tensor curAxialContactForceDerivative =
                                    curLineContacts[segI][nbI]
                                   .axialContactForceDerivative();

                                const vector curCircumContactForce =
                                    curLineContacts[segI][nbI]
                                   .circumferentialContactForce();

                                const tensor curCircumContactForceDerivative =
                                    curLineContacts[segI][nbI]
                                   .circumferentialContactForceDerivative();

                                const tensor curCircumContactForceThetaDerivative =
                                    curLineContacts[segI][nbI]
                                   .circumferentialContactForceThetaDerivative();

                                const tensor curContactMomentDerivative =
                                    curLineContacts[segI][nbI]
                                   .circumferentialContactMomentDerivative();

                                const vector dRdZeta =
                                    contact().splines()[nbI]
                                   .paramFirstDerivative(neiSegI, neiSegParam);
                                   
                                // What is happening here? Check with ZT   

                                tensor ImNN = tensor::zero;
                                vector N = curContactForce;
                                if (mag(N) > SMALL)
                                {
                                    N /= mag(N);
                                    ImNN = tensor::I - (N*N);
                                }

                                const tensor TT = //tensor::zero;
                                    (dRdZeta*dRdZeta)/((dRdZeta & dRdZeta));

                                curContactForceDerivative +=
                                (
                                    mag(curContactForce)*(ImNN - TT)
                                   *L()[cellID]/curContactDistance
                                );

                                // Add frictional component
                                if (mag(curCircumContactForce) > SMALL)
                                {
                                    curContactForceDerivative +=
                                        curCircumContactForceDerivative
                                       *L()[cellID];

                                    curContactForceDerivative +=
                                        curAxialContactForceDerivative
                                       *L()[cellID];
                                }

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

                                //- Contact force derivative (Theta)
                                diag[cellID](0,3) +=
                                    curCircumContactForceThetaDerivative.xx()
                                   *L()[cellID];
                                diag[cellID](0,4) +=
                                    curCircumContactForceThetaDerivative.xy()
                                   *L()[cellID];
                                diag[cellID](0,5) +=
                                    curCircumContactForceThetaDerivative.xz()
                                   *L()[cellID];

                                diag[cellID](1,3) +=
                                    curCircumContactForceThetaDerivative.yx()
                                   *L()[cellID];
                                diag[cellID](1,4) +=
                                    curCircumContactForceThetaDerivative.yy()
                                   *L()[cellID];
                                diag[cellID](1,5) +=
                                    curCircumContactForceThetaDerivative.yz()
                                   *L()[cellID];

                                diag[cellID](2,3) +=
                                    curCircumContactForceThetaDerivative.zx()
                                   *L()[cellID];
                                diag[cellID](2,4) +=
                                    curCircumContactForceThetaDerivative.zy()
                                   *L()[cellID];
                                diag[cellID](2,5) +=
                                    curCircumContactForceThetaDerivative.zz()
                                   *L()[cellID];
                            
                                //- Contact moment derivative
                                diag[cellID](3,3) +=
                                    curContactMomentDerivative.xx()
                                   *L()[cellID];
                                diag[cellID](3,4) +=
                                    curContactMomentDerivative.xy()
                                   *L()[cellID];
                                diag[cellID](3,5) +=
                                    curContactMomentDerivative.xz()
                                   *L()[cellID];
                                
                                diag[cellID](4,3) +=
                                    curContactMomentDerivative.yx()
                                   *L()[cellID];
                                diag[cellID](4,4) +=
                                    curContactMomentDerivative.yy()
                                   *L()[cellID];
                                diag[cellID](4,5) +=
                                    curContactMomentDerivative.yz()
                                   *L()[cellID];

                                diag[cellID](5,3) +=
                                    curContactMomentDerivative.zx()
                                   *L()[cellID];
                                diag[cellID](5,4) +=
                                    curContactMomentDerivative.zy()
                                   *L()[cellID];
                                diag[cellID](5,5) +=
                                    curContactMomentDerivative.zz()
                                   *L()[cellID];

                                // Off-diagonal

                                scalar w0 =
                                    curLineContacts[segI][nbI].secondBeamWeight();
                                scalar w1 = 1.0 - w0;

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


                                // Contact force derivative (Theta)

                                //- lower
                                lcCoeffs[bI][segI][nbI].first()(0, 3) +=
                                    curCircumContactForceThetaDerivative.xx()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(0, 4) +=
                                    curCircumContactForceThetaDerivative.xy()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(0, 5) +=
                                    curCircumContactForceThetaDerivative.xz()
                                   *L()[cellID]*w0;

                                lcCoeffs[bI][segI][nbI].first()(1, 3) +=
                                    curCircumContactForceThetaDerivative.yx()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(1, 4) +=
                                    curCircumContactForceThetaDerivative.yy()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(1, 5) +=
                                    curCircumContactForceThetaDerivative.yz()
                                   *L()[cellID]*w0;

                                lcCoeffs[bI][segI][nbI].first()(2, 3) +=
                                    curCircumContactForceThetaDerivative.zx()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(2, 4) +=
                                    curCircumContactForceThetaDerivative.zy()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(2, 5) +=
                                    curCircumContactForceThetaDerivative.zz()
                                   *L()[cellID]*w0;

                                //- upper
                                lcCoeffs[bI][segI][nbI].second()(0, 3) +=
                                    curCircumContactForceThetaDerivative.xx()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(0, 4) +=
                                    curCircumContactForceThetaDerivative.xy()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(0, 5) +=
                                    curCircumContactForceThetaDerivative.xz()
                                   *L()[cellID]*w1;

                                lcCoeffs[bI][segI][nbI].second()(1, 3) +=
                                    curCircumContactForceThetaDerivative.yx()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(1, 4) +=
                                    curCircumContactForceThetaDerivative.yy()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(1, 5) +=
                                    curCircumContactForceThetaDerivative.yz()
                                   *L()[cellID]*w1;

                                lcCoeffs[bI][segI][nbI].second()(2, 3) +=
                                    curCircumContactForceThetaDerivative.zx()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(2, 4) +=
                                    curCircumContactForceThetaDerivative.zy()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(2, 5) +=
                                    curCircumContactForceThetaDerivative.zz()
                                   *L()[cellID]*w1;

                                // Contact moment derivative

                                //- lower
                                lcCoeffs[bI][segI][nbI].first()(3, 3) +=
                                    curContactMomentDerivative.xx()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(3, 4) +=
                                    curContactMomentDerivative.xy()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(3, 5) +=
                                    curContactMomentDerivative.xz()
                                   *L()[cellID]*w0;

                                lcCoeffs[bI][segI][nbI].first()(4, 3) +=
                                    curContactMomentDerivative.yx()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(4, 4) +=
                                    curContactMomentDerivative.yy()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(4, 5) +=
                                    curContactMomentDerivative.yz()
                                   *L()[cellID]*w0;

                                lcCoeffs[bI][segI][nbI].first()(5, 3) +=
                                    curContactMomentDerivative.zx()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(5, 4) +=
                                    curContactMomentDerivative.zy()
                                   *L()[cellID]*w0;
                                lcCoeffs[bI][segI][nbI].first()(5, 5) +=
                                    curContactMomentDerivative.zz()
                                   *L()[cellID]*w0;

                                //- upper
                                lcCoeffs[bI][segI][nbI].second()(3, 3) +=
                                    curContactMomentDerivative.xx()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(3, 4) +=
                                    curContactMomentDerivative.xy()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(3, 5) +=
                                    curContactMomentDerivative.xz()
                                   *L()[cellID]*w1;

                                lcCoeffs[bI][segI][nbI].second()(4, 3) +=
                                    curContactMomentDerivative.yx()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(4, 4) +=
                                    curContactMomentDerivative.yy()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(4, 5) +=
                                    curContactMomentDerivative.yz()
                                   *L()[cellID]*w1;

                                lcCoeffs[bI][segI][nbI].second()(5, 3) +=
                                    curContactMomentDerivative.zx()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(5, 4) +=
                                    curContactMomentDerivative.zy()
                                   *L()[cellID]*w1;
                                lcCoeffs[bI][segI][nbI].second()(5, 5) +=
                                    curContactMomentDerivative.zz()
                                   *L()[cellID]*w1;
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
    multibeamFvBlockMatrix& eqn
)
{
    if
    (
        (mesh().cellZones().size() > 1)
     && !Pstream::parRun()
    )
    {
    Info << "Inside applyPointContact() in applyContacts.C \n" << endl;
    // Get matrix diagonal
    tensor6Field& diag = eqn.diag().asSquare();
        
    // Get matrix upper coeffs
    tensor6Field& upper = eqn.upper().asSquare();
    
    // Get matrix lower coeffs
    tensor6Field& lower = eqn.lower().asSquare();
    
    // Grab source
    vector6Field& source = eqn.source();

    // Grab beam coupling coeffs
    tensor6PairListList& pcCoeffs = eqn.pointContactCoeffs();
    
    // Info << "pcCoeffs: " << pcCoeffs << endl;

    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    const surfaceVectorField dRdS = dR0Ds_ + fvc::snGrad(W_);    

    const surfaceScalarField& dc = mesh().deltaCoeffs();

    // label start = 0;
    forAll(contact().pointContacts(), pcI)
    {
	Info << "Inside contact().pointContacts() in applyContacts.C" << endl;
        labelList bI(2, -1);
        labelList globalSegI(2, -1);
        scalarField zeta(2, -2);
        vectorList contactForce(2, vector::zero);
        vectorList tangContactForce(2, vector::zero);

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
            
            contactForce[0] =
                contact().pointContacts()[pcI].normalContactForce()
              + contact().pointContacts()[pcI].firstBeamTangContactForce()
              - contact().pointContacts()[pcI].secondBeamTangContactForce();
	    /*
	    Info << "\n bI[0]  " << bI[0] 
		 << "\n segI " << segI 
		 << "\n zeta[0] " << zeta[0]
		 << "\n globalSegI[0] " << globalSegI[0]
		 << "\n globalSegI[1] " << globalSegI[1]
		 << "\n normal con f " << contact().pointContacts()[pcI].normalContactForce()
		 << "\n fir B tan f " << contact().pointContacts()[pcI].firstBeamTangContactForce()
		 << "\n sec B tan f " << contact().pointContacts()[pcI].secondBeamTangContactForce()
		 << endl;
		 */
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
                
            contactForce[1] = -contactForce[0];
        }

        forAll(globalSegI, sI)
        {
            vector DR = vector::zero;

            if (zeta[sI] > 0)
            {
                label faceID = findIndex(own, globalSegI[sI]);
                if (faceID == -1) // last cell
                {
                    const unallocLabelList& faceCells =
                        mesh().boundary()[endPatchIndex(bI[sI])].faceCells();

                    label bFaceID = findIndex(faceCells, globalSegI[sI]);
                
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
                label faceID = findIndex(nei, globalSegI[sI]);
                if (faceID == -1) // first cell
                {
                    const unallocLabelList& faceCells =
                        mesh().boundary()[startPatchIndex(bI[sI])].faceCells();

                    label bFaceID = findIndex(faceCells, globalSegI[sI]);

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

            vector F0 = contactForce[sI];
	    
	    Info << "Contact Force F0 explicit term in applyContacts.C file "
		 << tab << F0 << endl;
		 
            vector M0 =  (spinTensor(DR) & F0);
	    //Info << "Contact  moment: " << M0 << endl;

            // label index = 6*globalSegI[sI];
            // b(index++) -= F0.x();
            // b(index++) -= F0.y();
            // b(index++) -= F0.z();
            // b(index++) -= M0.x();
            // b(index++) -= M0.y();
            // b(index++) -= M0.z();

            // W equation
            source[globalSegI[sI]](0) -= F0.x();
            source[globalSegI[sI]](1) -= F0.y();
            source[globalSegI[sI]](2) -= F0.z();

            // Theta equation
            source[globalSegI[sI]](3) -= M0.x();
            source[globalSegI[sI]](4) -= M0.y();
            source[globalSegI[sI]](5) -= M0.z();
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
	// SB
	labelList lenZetaI(2, -1);
	labelList lenNeiZetaI(2, -1);
	//--//
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
        vectorList contactForce(2, vector::zero);
        tensorList contactForceDerivative(2, tensor::zero);

        vectorList tangentialContactForce(2, vector::zero);
        tensorList tangContactForceDerivative(2, tensor::zero);
	

        //Owner
        {
            bI[0] = contact().pointContacts()[pcI].firstBeam();
            segI[0] = contact().pointContacts()[pcI].firstBeamSegment();
            zeta[0] = contact().pointContacts()[pcI].firstBeamZeta();
	    // SB
	    //scalarField L = contact().splines()[bI[segI[0]]].segLengths();
	    //Info << "bI " << bI[0] << "segI " << segI[0] << "zeta " << zeta[0] << endl;
	    //lenZetaI[0] = 0.5*mag(zeta[0])*L[segI[0]];
	    //Info << "abs(zeta) " << mag(zeta[0]) << endl;
	    //Info << "len[segI] " << L[segI[0]] << endl;
	    //Info << "LenZeta " << lenZetaI[0] << endl;
	    
	    //--//
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
            contactForce[0] =
                contact().pointContacts()[pcI].normalContactForce();
            contactForceDerivative[0] =
                contact().pointContacts()[pcI].normalContactForceDerivative();

            tangentialContactForce[0] = 
                contact().pointContacts()[pcI].firstBeamTangContactForce();

            tangContactForceDerivative[0] = 
                contact().pointContacts()[pcI]
               .firstBeamTangContactForceDerivative()
              - contact().pointContacts()[pcI]
               .secondBeamTangContactForceDerivative();

	    Info << "contact distance: " 
	      << contact().pointContacts()[pcI].delta()
	      << endl;
	 
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
            contactForce[1] = -contactForce[0];
	    
            contactForceDerivative[1] = contactForceDerivative[0];

            tangentialContactForce[1] = 
                contact().pointContacts()[pcI]
               .secondBeamTangContactForce();
                      
            tangContactForceDerivative[1] = 
                tangContactForceDerivative[0];
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
            vector DR = vector::zero;
            if (zeta[sI] > 0)
            {
                label faceID = findIndex(own, globalSegI[sI]);
                if (faceID == -1) // last cell
                {
                    const unallocLabelList& faceCells =
                        mesh().boundary()[endPatchIndex(bI[sI])].faceCells();

                    label bFaceID =
                        findIndex(faceCells, globalSegI[sI]);

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
                label faceID = findIndex(nei, globalSegI[sI]);
                if (faceID == -1) // first cell
                {
                    const unallocLabelList& faceCells =
                        mesh().boundary()[startPatchIndex(bI[sI])].faceCells();

                    label bFaceID =
                        findIndex(faceCells, globalSegI[sI]);

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

            vector N = contactForce[sI];
            tensor ImNN = tensor::zero;
            if (mag(N) > SMALL)
            {
                N /= mag(N);
                ImNN = tensor::I - (N*N);
            }

            //
            vector dRdZeta =
                contact().splines()[bI[sI]]
               .paramFirstDerivative(segI[sI], zeta[sI]);
            vector neiDRdZeta =
                contact().splines()[nbI[sI]]
               .paramFirstDerivative(neiSegI[sI], neiZeta[sI]);
            vector dR2dZeta2 =
                contact().splines()[bI[sI]]
               .paramSecondDerivative(segI[sI], zeta[sI]);
            vector neiDR2dZeta2 =
                contact().splines()[nbI[sI]]
               .paramSecondDerivative(neiSegI[sI], neiZeta[sI]);

            vector delta = contactDistance[sI]*N;

            scalarSquareMatrix M(2, 0.0);
	    
	    //Original expressions
	     
            M[0][0] = (dRdZeta & neiDRdZeta);
            M[0][1] = ((delta & neiDR2dZeta2) - (neiDRdZeta & neiDRdZeta));
            M[1][0] = ((delta & dR2dZeta2) + (dRdZeta & dRdZeta));
            M[1][1] = -(neiDRdZeta & dRdZeta);

	    
            scalarSquareMatrix invM = M.LUinvert();

            tensor TT = 
              - (dRdZeta*(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta))
              + (neiDRdZeta*(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta));
	    

            tensor FcDn =
            (
                mag(contactForce[sI])
               *(ImNN + (ImNN & TT))/contactDistance[sI]
            );
	    
	    
	    // Implicit part of moment is turned off because it does not
	    // make the convergence better
            tensor DMCoeff =  tensor::zero;
       //     (
       //         spinTensor(DR)
       //      & (contactForceDerivative[sI] + FcDn)
       //     );
 
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
              + wDiag*FcDn.xx();
            diag[diagCell](0, 1) +=
                wDiag*contactForceDerivative[sI].xy()
              + wDiag*tangContactForceDerivative[sI].xy()
              + wDiag*FcDn.xy();
            diag[diagCell](0, 2) +=
                wDiag*contactForceDerivative[sI].xz()
              + wDiag*tangContactForceDerivative[sI].xz()
              + wDiag*FcDn.xz();

            diag[diagCell](1, 0) +=
                wDiag*contactForceDerivative[sI].yx()
              + wDiag*tangContactForceDerivative[sI].yx()
              + wDiag*FcDn.yx();      
            diag[diagCell](1, 1) +=
                wDiag*contactForceDerivative[sI].yy()
              + wDiag*tangContactForceDerivative[sI].yy()
              + wDiag*FcDn.yy();      
            diag[diagCell](1, 2) +=
                wDiag*contactForceDerivative[sI].yz()
              + wDiag*tangContactForceDerivative[sI].yz()
              + wDiag*FcDn.yz();

            diag[diagCell](2, 0) +=
                wDiag*contactForceDerivative[sI].zx()
              + wDiag*tangContactForceDerivative[sI].zx()
              + wDiag*FcDn.zx();
            diag[diagCell](2, 1) +=
                wDiag*contactForceDerivative[sI].zy()
              + wDiag*tangContactForceDerivative[sI].zy()
              + wDiag*FcDn.zy();
            diag[diagCell](2, 2) +=
                wDiag*contactForceDerivative[sI].zz()
              + wDiag*tangContactForceDerivative[sI].zz()
              + wDiag*FcDn.zz();

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
                tensor6* curCoeffPtr = &(upper[sharedFaceIndex]);
                if (own[sharedFaceIndex] != globalSegI[sI])
                {
                    curCoeffPtr = &(lower[sharedFaceIndex]);
                }
                tensor6& curCoeff = *curCoeffPtr;
		
                //- Force
                curCoeff(0, 0) +=
                    (1.0-wDiag)*contactForceDerivative[sI].xx()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].xx()
                  + (1.0-wDiag)*FcDn.xx();
                curCoeff(0, 1) +=
                    (1.0-wDiag)*contactForceDerivative[sI].xy()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].xy()
                  + (1.0-wDiag)*FcDn.xy();
                curCoeff(0, 2) +=
                    (1.0-wDiag)*contactForceDerivative[sI].xz()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].xz()
                  + (1.0-wDiag)*FcDn.xz();

                curCoeff(1, 0) +=
                    (1.0-wDiag)*contactForceDerivative[sI].yx()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].yx()
                  + (1.0-wDiag)*FcDn.yx();
                curCoeff(1, 1) +=
                    (1.0-wDiag)*contactForceDerivative[sI].yy()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].yy()
                  + (1.0-wDiag)*FcDn.yy();
                curCoeff(1, 2) +=
                    (1.0-wDiag)*contactForceDerivative[sI].yz()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].yz()
                  + (1.0-wDiag)*FcDn.yz();

                curCoeff(2, 0) +=
                    (1.0-wDiag)*contactForceDerivative[sI].zx()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].zx()
                  + (1.0-wDiag)*FcDn.zx();
                curCoeff(2, 1) +=
                    (1.0-wDiag)*contactForceDerivative[sI].zy()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].zy()
                  + (1.0-wDiag)*FcDn.zy();
                curCoeff(2, 2) +=
                    (1.0-wDiag)*contactForceDerivative[sI].zz()
                  + (1.0-wDiag)*tangContactForceDerivative[sI].zz()
                  + (1.0-wDiag)*FcDn.zz();
                
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

            //-- IN0
            //- Force
            pcCoeffs[pcI][sI].first()(0,0) -=
                w0*contactForceDerivative[sI].xx()
              + w0*tangContactForceDerivative[sI].xx()
              + w0*FcDn.xx();
            pcCoeffs[pcI][sI].first()(0,1) -=
                w0*contactForceDerivative[sI].xy()
              + w0*tangContactForceDerivative[sI].xy()
              + w0*FcDn.xy();
            pcCoeffs[pcI][sI].first()(0,2) -=
                w0*contactForceDerivative[sI].xz()
              + w0*tangContactForceDerivative[sI].xz()
              + w0*FcDn.xz();

            pcCoeffs[pcI][sI].first()(1,0) -=
                w0*contactForceDerivative[sI].yx()
              + w0*tangContactForceDerivative[sI].yx()
              + w0*FcDn.yx();
            pcCoeffs[pcI][sI].first()(1,1) -=
                w0*contactForceDerivative[sI].yy()
              + w0*tangContactForceDerivative[sI].yy()
              + w0*FcDn.yy();
            pcCoeffs[pcI][sI].first()(1,2) -=
                w0*contactForceDerivative[sI].yz()
              + w0*tangContactForceDerivative[sI].yz()
              + w0*FcDn.yz();

            pcCoeffs[pcI][sI].first()(2,0) -=
                w0*contactForceDerivative[sI].zx()
              + w0*tangContactForceDerivative[sI].zx()
              + w0*FcDn.zx();
            pcCoeffs[pcI][sI].first()(2,1) -=
                w0*contactForceDerivative[sI].zy()
              + w0*tangContactForceDerivative[sI].zy()
              + w0*FcDn.zy();
            pcCoeffs[pcI][sI].first()(2,2) -=
                w0*contactForceDerivative[sI].zz()
              + w0*tangContactForceDerivative[sI].zz()
              + w0*FcDn.zz();
	    
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
              + w1*tangContactForceDerivative[sI].xx()
              + w1*FcDn.xx();
            pcCoeffs[pcI][sI].second()(0,1) -=
                w1*contactForceDerivative[sI].xy()
              + w1*tangContactForceDerivative[sI].xy()
              + w1*FcDn.xy();
            pcCoeffs[pcI][sI].second()(0,2) -=
                w1*contactForceDerivative[sI].xz()
              + w1*tangContactForceDerivative[sI].xz()
              + w1*FcDn.xz();

            pcCoeffs[pcI][sI].second()(1,0) -=
                w1*contactForceDerivative[sI].yx()
              + w1*tangContactForceDerivative[sI].yx()
              + w1*FcDn.yx();
            pcCoeffs[pcI][sI].second()(1,1) -=
                w1*contactForceDerivative[sI].yy()
              + w1*tangContactForceDerivative[sI].yy()
              + w1*FcDn.yy();
            pcCoeffs[pcI][sI].second()(1,2) -=
                w1*contactForceDerivative[sI].yz()
              + w1*tangContactForceDerivative[sI].yz()
              + w1*FcDn.yz();

            pcCoeffs[pcI][sI].second()(2,0) -=
                w1*contactForceDerivative[sI].zx()
              + w1*tangContactForceDerivative[sI].zx()
              + w1*FcDn.zx();
            pcCoeffs[pcI][sI].second()(2,1) -=
                w1*contactForceDerivative[sI].zy()
              + w1*tangContactForceDerivative[sI].zy()
              + w1*FcDn.zy();
            pcCoeffs[pcI][sI].second()(2,2) -=
                w1*contactForceDerivative[sI].zz()
              + w1*tangContactForceDerivative[sI].zz()
              + w1*FcDn.zz();

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
            pcCoeffs[pcI][sI].second()(5,0) -=
                w1*DMCoeff.zy();
            pcCoeffs[pcI][sI].second()(5,0) -=
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
    
void coupledTotalLagNewtonRaphsonBeam::applyConicalPulleysContact
(
    multibeamFvBlockMatrix& eqn
)
{
    if (conicalPulleys().size())
    {
        // Get matrix diagonal
        tensor6Field& blockDiag = eqn.diag().asSquare();

        // Grab source
        vector6Field& source = eqn.source();

        // Grab off-diagonal 
        tensor6Field& upper = eqn.upper().asSquare();

        // Grab off-diagonal
        tensor6Field& lower = eqn.lower().asSquare();
    
        // Info << "Add explicit contribution of conical pulley contac force"
        //      << endl;

        label nPulleys = conicalPulleys().size();

        label start = 0;
        // label index = 0;
    
        label nBeams = contact().splines().size();
        // const scalarField& L = L();

        // const cellList& cells = mesh().cells();
        // const labelListList& cc = mesh().cellCells();

        // const volScalarField& acpc =
        //     contact().activeConicalPulleyContacts();

        for (label bI=0; bI<nBeams; bI++)
        {
            const conicalPulleyContactListList& curPulleyContacts =
                contact().conicalPulleyContacts()[bI];

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
            
            for
            (
                label segI=0;
                segI<contact().splines()[bI].nSegments();
                segI++
            )
            {
                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    for (label cI=0; cI<2; cI++)
                    {
                        curq[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .normalContactForce()[cI];

                        curq[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .frictionalContactForce()[cI];

                        // Contact force due to friction in pulley bearings
                        curq[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .axialBearingFrictionForce()[cI];
                    
                        curm[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .frictionalContactMoment()[cI];
                    }
                }
            }

            forAll(curq, segI)
            {
                label globalSegI = start + segI;

                // W equation
                source[globalSegI](0) -= curq[segI].x()*L()[globalSegI];
                source[globalSegI](1) -= curq[segI].y()*L()[globalSegI];
                source[globalSegI](2) -= curq[segI].z()*L()[globalSegI];

                // Theta equation
                source[globalSegI](3) -= curm[segI].x()*L()[globalSegI];
                source[globalSegI](4) -= curm[segI].y()*L()[globalSegI];
                source[globalSegI](5) -= curm[segI].z()*L()[globalSegI];

                // b(index++) -= curm[segI].x()*L[start+segI];
                // b(index++) -= curm[segI].y()*L[start+segI];
                // b(index++) -= curm[segI].z()*L[start+segI];
                // index++;
                // index++;
                // index++;
            }

            start += curq.size();
        }

        // Info << "Add implicit contribution of conical pulley contac force"
        //      << endl;
        
        const labelList& upperNeiFaces = upperNeiCellFaces();
        const labelList& lowerNeiFaces = lowerNeiCellFaces();
            
        start = 0;
        for (label bI=0; bI<nBeams; bI++)
        {
            const conicalPulleyContactListList& curPulleyContacts =
                contact().conicalPulleyContacts()[bI];
    
            // const vectorField& dRdS = contact().splines()[bI].dRdS();

            label nSeg = contact().splines()[bI].nSegments();

            for (label segI=0; segI<nSeg; segI++)
            {
                label globalSegI = start + segI;

                vector T = contact().splines()[bI].dRdS()[segI];
                T /= mag(T) + SMALL;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    const tensorField& curContactForceDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForceDerivative();

                    const tensorField& curContactForceDirectionDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForceDirectionDerivative();

                    tensorField FcDn =
                        curContactForceDirectionDerivative
                       *L()[globalSegI];
                    
                    const vectorField& curContactMoment =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactMoment();

                    const tensorField& curContactMomentDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactMomentDerivative();

                    const tensorField& curContactMomentDerivativeOverFcn =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactMomentDerivativeOverFcn();

                    const tensorField& curFrictionalContactForceDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactForceDerivative();

                    const tensorField& curFrictionalContactForceDerivativeOverFcn =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactForceDerivativeOverFcn();

                    for (label cI=0; cI<2; cI++)
                    {
                        //- Contact force derivative (W)
                        blockDiag[globalSegI](0,0) +=
                            curContactForceDerivative[cI].xx()*L()[globalSegI]
                          + FcDn[cI].xx();
                        blockDiag[globalSegI](0,1) +=
                            curContactForceDerivative[cI].xy()*L()[globalSegI]
                          + FcDn[cI].xy();
                        blockDiag[globalSegI](0,2) +=
                            curContactForceDerivative[cI].xz()*L()[globalSegI]
                          + FcDn[cI].xz();

                        blockDiag[globalSegI](1,0) +=
                            curContactForceDerivative[cI].yx()*L()[globalSegI]
                          + FcDn[cI].yx();
                        blockDiag[globalSegI](1,1) +=
                            curContactForceDerivative[cI].yy()*L()[globalSegI]
                          + FcDn[cI].yy();
                        blockDiag[globalSegI](1,2) +=
                            curContactForceDerivative[cI].yz()*L()[globalSegI]
                          + FcDn[cI].yz();

                        blockDiag[globalSegI](2,0) +=
                            curContactForceDerivative[cI].zx()*L()[globalSegI]
                          + FcDn[cI].zx();
                        blockDiag[globalSegI](2,1) +=
                            curContactForceDerivative[cI].zy()*L()[globalSegI]
                          + FcDn[cI].zy();
                        blockDiag[globalSegI](2,2) +=
                            curContactForceDerivative[cI].zz()*L()[globalSegI]
                          + FcDn[cI].zz();

                        //- Frictional contact moment derivative (Theta)
                        blockDiag[globalSegI](3,3) +=
                            curContactMomentDerivative[cI].xx()*L()[globalSegI];
                        blockDiag[globalSegI](3,4) +=
                            curContactMomentDerivative[cI].xy()*L()[globalSegI];
                        blockDiag[globalSegI](3,5) +=
                            curContactMomentDerivative[cI].xz()*L()[globalSegI];

                        blockDiag[globalSegI](4,3) +=
                            curContactMomentDerivative[cI].yx()*L()[globalSegI];
                        blockDiag[globalSegI](4,4) +=
                            curContactMomentDerivative[cI].yy()*L()[globalSegI];
                        blockDiag[globalSegI](4,5) +=
                            curContactMomentDerivative[cI].yz()*L()[globalSegI];

                        blockDiag[globalSegI](5,3) +=
                            curContactMomentDerivative[cI].zx()*L()[globalSegI];
                        blockDiag[globalSegI](5,4) +=
                            curContactMomentDerivative[cI].zy()*L()[globalSegI];
                        blockDiag[globalSegI](5,5) +=
                            curContactMomentDerivative[cI].zz()*L()[globalSegI];

                        //- Frictional contact moment derivative (W)
                        blockDiag[globalSegI](3,0) +=
                            curContactMomentDerivativeOverFcn[cI].xx()*L()[globalSegI];
                        blockDiag[globalSegI](3,1) +=
                            curContactMomentDerivativeOverFcn[cI].xy()*L()[globalSegI];
                        blockDiag[globalSegI](3,2) +=
                            curContactMomentDerivativeOverFcn[cI].xz()*L()[globalSegI];

                        blockDiag[globalSegI](4,0) +=
                            curContactMomentDerivativeOverFcn[cI].yx()*L()[globalSegI];
                        blockDiag[globalSegI](4,1) +=
                            curContactMomentDerivativeOverFcn[cI].yy()*L()[globalSegI];
                        blockDiag[globalSegI](4,2) +=
                            curContactMomentDerivativeOverFcn[cI].yz()*L()[globalSegI];

                        blockDiag[globalSegI](5,0) +=
                            curContactMomentDerivativeOverFcn[cI].zx()*L()[globalSegI];
                        blockDiag[globalSegI](5,1) +=
                            curContactMomentDerivativeOverFcn[cI].zy()*L()[globalSegI];
                        blockDiag[globalSegI](5,2) +=
                            curContactMomentDerivativeOverFcn[cI].zz()*L()[globalSegI];

                        //- Frictional contact force derivative (Theta)
                        if (true)
                        {
                        blockDiag[globalSegI](0,3) +=
                            curFrictionalContactForceDerivative[cI].xx()*L()[globalSegI];
                        blockDiag[globalSegI](0,4) +=
                            curFrictionalContactForceDerivative[cI].xy()*L()[globalSegI];
                        blockDiag[globalSegI](0,5) +=
                            curFrictionalContactForceDerivative[cI].xz()*L()[globalSegI];

                        blockDiag[globalSegI](1,3) +=
                            curFrictionalContactForceDerivative[cI].yx()*L()[globalSegI];
                        blockDiag[globalSegI](1,4) +=
                            curFrictionalContactForceDerivative[cI].yy()*L()[globalSegI];
                        blockDiag[globalSegI](1,5) +=
                            curFrictionalContactForceDerivative[cI].yz()*L()[globalSegI];

                        blockDiag[globalSegI](2,3) +=
                            curFrictionalContactForceDerivative[cI].zx()*L()[globalSegI];
                        blockDiag[globalSegI](2,4) +=
                            curFrictionalContactForceDerivative[cI].zy()*L()[globalSegI];
                        blockDiag[globalSegI](2,5) +=
                            curFrictionalContactForceDerivative[cI].zz()*L()[globalSegI];
                        }

                        //- Frictional contact force derivative (W)
                        if (true)
                        {
                        blockDiag[globalSegI](0,0) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].xx()*L()[globalSegI];
                        blockDiag[globalSegI](0,1) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].xy()*L()[globalSegI];
                        blockDiag[globalSegI](0,2) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].xz()*L()[globalSegI];

                        blockDiag[globalSegI](1,0) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].yx()*L()[globalSegI];
                        blockDiag[globalSegI](1,1) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].yy()*L()[globalSegI];
                        blockDiag[globalSegI](1,2) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].yz()*L()[globalSegI];

                        blockDiag[globalSegI](2,0) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].zx()*L()[globalSegI];
                        blockDiag[globalSegI](2,1) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].zy()*L()[globalSegI];
                        blockDiag[globalSegI](2,2) +=
                            curFrictionalContactForceDerivativeOverFcn[cI].zz()*L()[globalSegI];
                        }
                        
                        //- Contact moment direction derivative
                        if (false)
                        // if (cc[globalSegI].size() == 2)
                        {
                            const label upperNeiFace =
                                upperNeiFaces[globalSegI];
                            const label lowerNeiFace =
                                lowerNeiFaces[globalSegI];
                            
                            if (upperNeiFace < 0 || lowerNeiFace < 0)
                            {
                                FatalErrorIn
                                (
                                    "coupledTotalLagNewtonRaphsonBeam::"
                                    "applyConicalPulleysContact(...)"
                                )
                                    << "Can not find upper and lower neighbour."
                                    << abort(FatalError);
                            }

                            tensor curCoeff = 
                                tensor::I*(curContactMoment[cI] & T)/2;
                            // already multiplied by L()[globalSegI]

                            upper[upperNeiFace](3,0) += curCoeff.xx();
                            upper[upperNeiFace](3,1) += curCoeff.xy();
                            upper[upperNeiFace](3,2) += curCoeff.xz();

                            upper[upperNeiFace](4,0) += curCoeff.yx();
                            upper[upperNeiFace](4,1) += curCoeff.yy();
                            upper[upperNeiFace](4,2) += curCoeff.yz();

                            upper[upperNeiFace](5,0) += curCoeff.zx();
                            upper[upperNeiFace](5,1) += curCoeff.zy();
                            upper[upperNeiFace](5,2) += curCoeff.zz();

                            lower[lowerNeiFace](3,0) -= curCoeff.xx();
                            lower[lowerNeiFace](3,1) -= curCoeff.xy();
                            lower[lowerNeiFace](3,2) -= curCoeff.xz();

                            lower[lowerNeiFace](4,0) -= curCoeff.yx();
                            lower[lowerNeiFace](4,1) -= curCoeff.yy();
                            lower[lowerNeiFace](4,2) -= curCoeff.yz();

                            lower[lowerNeiFace](5,0) -= curCoeff.zx();
                            lower[lowerNeiFace](5,1) -= curCoeff.zy();
                            lower[lowerNeiFace](5,2) -= curCoeff.zz();
                        }

                        if (false)
                        {
                            tensor curCoeff =
                               -(curContactMoment[cI] & T) * spinTensor(T);

                            blockDiag[globalSegI](3,3) +=
                                curCoeff.xx()*L()[globalSegI];
                            blockDiag[globalSegI](3,4) +=
                                curCoeff.xy()*L()[globalSegI];
                            blockDiag[globalSegI](3,5) +=
                                curCoeff.xz()*L()[globalSegI];

                            blockDiag[globalSegI](4,3) +=
                                curCoeff.yx()*L()[globalSegI];
                            blockDiag[globalSegI](4,4) +=
                                curCoeff.yy()*L()[globalSegI];
                            blockDiag[globalSegI](4,5) +=
                                curCoeff.yz()*L()[globalSegI];

                            blockDiag[globalSegI](5,3) +=
                                curCoeff.zx()*L()[globalSegI];
                            blockDiag[globalSegI](5,4) +=
                                curCoeff.zy()*L()[globalSegI];
                            blockDiag[globalSegI](5,5) +=
                                curCoeff.zz()*L()[globalSegI];
                        }
                    }
                }
            }

            start += nSeg;
        }

        // if (isA<lagrangeMultipliers>(contact().normalModel()))
        // {
        //     #include "applyLagrangeMultiplierConicalPulleyContact.H"
        // }
    }
}



void coupledTotalLagNewtonRaphsonBeam::applyConicalPulleysContactNew
(
    multibeamFvBlockMatrix& eqn
)
{
    if (conicalPulleys().size())
    {
        // Get matrix diagonal
        tensor6Field& blockDiag = eqn.diag().asSquare();

        // Grab source
        vector6Field& source = eqn.source();

        // Grab off-diagonal 
        tensor6Field& upper = eqn.upper().asSquare();

        // Grab off-diagonal
        // tensor6Field& lower = eqn.lower().asSquare();
    
        // Info << "Add explicit contribution of conical pulley contac force"
        //      << endl;

        label nPulleys = conicalPulleys().size();

        label start = 0;
        // label index = 0;
    
        label nBeams = contact().splines().size();
        // const scalarField& L = L();

        // const cellList& cells = mesh().cells();
        // const labelListList& cc = mesh().cellCells();

        // const volScalarField& acpc =
        //     contact().activeConicalPulleyContacts();

        for (label bI=0; bI<nBeams; bI++)
        {
            const conicalPulleyContactListList& curPulleyContacts =
                contact().conicalPulleyContacts()[bI];

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
            
            for
            (
                label segI=0;
                segI<contact().splines()[bI].nSegments();
                segI++
            )
            {
                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    for (label cI=0; cI<2; cI++)
                    {
                        curq[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .normalContactForce()[cI];

                        curq[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .axialFrictionalContactForce()[cI];

                        curq[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .transversalFrictionalContactForce()[cI];

                        // // Contact force due to friction in pulley bearings
                        // curq[segI] +=
                        //     curPulleyContacts[segI][pulleyI]
                        //    .axialBearingFrictionForce()[cI];
                    
                        curm[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .transversalFrictionalContactMoment()[cI];
                    }
                }
            }

            forAll(curq, segI)
            {
                label globalSegI = start + segI;

                // W equation
                source[globalSegI](0) -= curq[segI].x()*L()[globalSegI];
                source[globalSegI](1) -= curq[segI].y()*L()[globalSegI];
                source[globalSegI](2) -= curq[segI].z()*L()[globalSegI];

                // Theta equation
                source[globalSegI](3) -= curm[segI].x()*L()[globalSegI];
                source[globalSegI](4) -= curm[segI].y()*L()[globalSegI];
                source[globalSegI](5) -= curm[segI].z()*L()[globalSegI];

                // b(index++) -= curm[segI].x()*L[start+segI];
                // b(index++) -= curm[segI].y()*L[start+segI];
                // b(index++) -= curm[segI].z()*L[start+segI];
                // index++;
                // index++;
                // index++;
            }

            start += curq.size();
        }

        // Info << "Add implicit contribution of conical pulley contac force"
        //      << endl;
        
        const labelList& upperNeiFaces = upperNeiCellFaces();
        const labelList& lowerNeiFaces = lowerNeiCellFaces();
        
        const scalarField& dc = mesh().deltaCoeffs().internalField();
        
        start = 0;
        for (label bI=0; bI<nBeams; bI++)
        {
            const conicalPulleyContactListList& curPulleyContacts =
                contact().conicalPulleyContacts()[bI];
    
            // const vectorField& dRdS = contact().splines()[bI].dRdS();

            label nSeg = contact().splines()[bI].nSegments();

            for (label segI=0; segI<nSeg; segI++)
            {
                label globalSegI = start + segI;

                vector T = contact().splines()[bI].dRdS()[segI];
                T /= mag(T) + SMALL;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    const tensorField& curNormalContactForceDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForceDerivative();

                    const tensorField&
                        curNormalContactForceDirectionDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForceDirectionDerivative();

                    // tensorField FcDn =
                    //     curContactForceDirectionDerivative*L()[globalSegI];

                    const vectorField& curAxialContactForce =
                        curPulleyContacts[segI][pulleyI]
                       .axialFrictionalContactForce();

                    // const vectorField& curContactMoment =
                    //     curPulleyContacts[segI][pulleyI]
                    //    .transversalFrictionalContactMoment();

                    const tensorField& curContactMomentDerivativePerDW =
                        curPulleyContacts[segI][pulleyI]
                       .transversalFrictionalContactMomentDerivativePerDW();

                    const tensorField& curContactMomentDerivativePerDTheta =
                        curPulleyContacts[segI][pulleyI]
                       .transversalFrictionalContactMomentDerivativePerDTheta();

                    const tensorField& curTransversalContactForceDerivativePerDW =
                        curPulleyContacts[segI][pulleyI]
                       .transversalFrictionalContactForceDerivativePerDW();

                    const tensorField&
                        curTransversalContactForceDerivativePerDTheta =
                        curPulleyContacts[segI][pulleyI]
                       .transversalFrictionalContactForceDerivativePerDTheta();

                    const tensorField& curAxialContactForceDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .axialFrictionalContactForceDerivative();

                    
                    // const tensorField& curContactMomentDerivativeOverFcn =
                    //     curPulleyContacts[segI][pulleyI]
                    //    .frictionalContactMomentDerivativeOverFcn();

                    // const tensorField& curFrictionalContactForceDerivative =
                    //     curPulleyContacts[segI][pulleyI]
                    //    .frictionalContactForceDerivative();

                    // const tensorField& curFrictionalContactForceDerivativeOverFcn =
                    //     curPulleyContacts[segI][pulleyI]
                    //    .frictionalContactForceDerivativeOverFcn();


                    tensorField curContactForceDerivativePerDW =
                        curNormalContactForceDerivative
                      + curNormalContactForceDirectionDerivative
                      + curAxialContactForceDerivative
                      + curTransversalContactForceDerivativePerDW;
                    curContactForceDerivativePerDW *= L()[globalSegI];

                    tensorField curContactForceDerivativePerDTheta =
                        curTransversalContactForceDerivativePerDTheta;
                    curContactForceDerivativePerDTheta *= L()[globalSegI];

                    // tensorField curContactMomentDerivativePerDW =
                    //     curTransversalContactMomentDerivativePerDW;
                    // curContactMomentDerivativePerDw *= L()[globalSegI];


                    for (label cI=0; cI<2; cI++)
                    {
                        //- Contact force derivative (W)
                        blockDiag[globalSegI](0,0) +=
                            curContactForceDerivativePerDW[cI].xx();
                        blockDiag[globalSegI](0,1) +=
                            curContactForceDerivativePerDW[cI].xy();
                        blockDiag[globalSegI](0,2) +=
                            curContactForceDerivativePerDW[cI].xz();

                        blockDiag[globalSegI](1,0) +=
                            curContactForceDerivativePerDW[cI].yx();
                        blockDiag[globalSegI](1,1) +=
                            curContactForceDerivativePerDW[cI].yy();
                        blockDiag[globalSegI](1,2) +=
                            curContactForceDerivativePerDW[cI].yz();

                        blockDiag[globalSegI](2,0) +=
                            curContactForceDerivativePerDW[cI].zx();
                        blockDiag[globalSegI](2,1) +=
                            curContactForceDerivativePerDW[cI].zy();
                        blockDiag[globalSegI](2,2) +=
                            curContactForceDerivativePerDW[cI].zz();

                        
                        //- Contact force derivative (Theta)
                        blockDiag[globalSegI](0,3) +=
                            curContactForceDerivativePerDTheta[cI].xx();
                        blockDiag[globalSegI](0,4) +=
                            curContactForceDerivativePerDTheta[cI].xy();
                        blockDiag[globalSegI](0,5) +=
                            curContactForceDerivativePerDTheta[cI].xz();

                        blockDiag[globalSegI](1,3) +=
                            curContactForceDerivativePerDTheta[cI].yx();
                        blockDiag[globalSegI](1,4) +=
                            curContactForceDerivativePerDTheta[cI].yy();
                        blockDiag[globalSegI](1,5) +=
                            curContactForceDerivativePerDTheta[cI].yz();

                        blockDiag[globalSegI](2,3) +=
                            curContactForceDerivativePerDTheta[cI].zx();
                        blockDiag[globalSegI](2,4) +=
                            curContactForceDerivativePerDTheta[cI].zy();
                        blockDiag[globalSegI](2,5) +=
                            curContactForceDerivativePerDTheta[cI].zz();


                        //- Frictional contact moment derivative (W)
                        blockDiag[globalSegI](3,0) +=
                            curContactMomentDerivativePerDW[cI].xx()*
                            L()[globalSegI];
                        blockDiag[globalSegI](3,1) +=
                            curContactMomentDerivativePerDW[cI].xy()*
                            L()[globalSegI];
                        blockDiag[globalSegI](3,2) +=
                            curContactMomentDerivativePerDW[cI].xz()*
                            L()[globalSegI];

                        blockDiag[globalSegI](4,0) +=
                            curContactMomentDerivativePerDW[cI].yx()*
                            L()[globalSegI];
                        blockDiag[globalSegI](4,1) +=
                            curContactMomentDerivativePerDW[cI].yy()*
                            L()[globalSegI];
                        blockDiag[globalSegI](4,2) +=
                            curContactMomentDerivativePerDW[cI].yz()*
                            L()[globalSegI];

                        blockDiag[globalSegI](5,0) +=
                            curContactMomentDerivativePerDW[cI].zx()*
                            L()[globalSegI];
                        blockDiag[globalSegI](5,1) +=
                            curContactMomentDerivativePerDW[cI].zy()*
                            L()[globalSegI];
                        blockDiag[globalSegI](5,2) +=
                            curContactMomentDerivativePerDW[cI].zz()*
                            L()[globalSegI];

                        
                        //- Frictional contact moment derivative (Theta)
                        blockDiag[globalSegI](3,3) +=
                            curContactMomentDerivativePerDTheta[cI].xx()*
                            L()[globalSegI];
                        blockDiag[globalSegI](3,4) +=
                            curContactMomentDerivativePerDTheta[cI].xy()*
                            L()[globalSegI];
                        blockDiag[globalSegI](3,5) +=
                            curContactMomentDerivativePerDTheta[cI].xz()*
                            L()[globalSegI];

                        blockDiag[globalSegI](4,3) +=
                            curContactMomentDerivativePerDTheta[cI].yx()*
                            L()[globalSegI];
                        blockDiag[globalSegI](4,4) +=
                            curContactMomentDerivativePerDTheta[cI].yy()*
                            L()[globalSegI];
                        blockDiag[globalSegI](4,5) +=
                            curContactMomentDerivativePerDTheta[cI].yz()*
                            L()[globalSegI];

                        blockDiag[globalSegI](5,3) +=
                            curContactMomentDerivativePerDTheta[cI].zx()*
                            L()[globalSegI];
                        blockDiag[globalSegI](5,4) +=
                            curContactMomentDerivativePerDTheta[cI].zy()*
                            L()[globalSegI];
                        blockDiag[globalSegI](5,5) +=
                            curContactMomentDerivativePerDTheta[cI].zz()*
                            L()[globalSegI];

                        
                        // //- Frictional contact force derivative (Theta)
                        // if (false)
                        // {
                        // blockDiag[globalSegI](0,3) +=
                        //     curFrictionalContactForceDerivative[cI].xx()*L()[globalSegI];
                        // blockDiag[globalSegI](0,4) +=
                        //     curFrictionalContactForceDerivative[cI].xy()*L()[globalSegI];
                        // blockDiag[globalSegI](0,5) +=
                        //     curFrictionalContactForceDerivative[cI].xz()*L()[globalSegI];

                        // blockDiag[globalSegI](1,3) +=
                        //     curFrictionalContactForceDerivative[cI].yx()*L()[globalSegI];
                        // blockDiag[globalSegI](1,4) +=
                        //     curFrictionalContactForceDerivative[cI].yy()*L()[globalSegI];
                        // blockDiag[globalSegI](1,5) +=
                        //     curFrictionalContactForceDerivative[cI].yz()*L()[globalSegI];

                        // blockDiag[globalSegI](2,3) +=
                        //     curFrictionalContactForceDerivative[cI].zx()*L()[globalSegI];
                        // blockDiag[globalSegI](2,4) +=
                        //     curFrictionalContactForceDerivative[cI].zy()*L()[globalSegI];
                        // blockDiag[globalSegI](2,5) +=
                        //     curFrictionalContactForceDerivative[cI].zz()*L()[globalSegI];
                        // }

                        // //- Frictional contact force derivative (W)
                        // if (false)
                        // {
                        // blockDiag[globalSegI](0,0) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].xx()*L()[globalSegI];
                        // blockDiag[globalSegI](0,1) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].xy()*L()[globalSegI];
                        // blockDiag[globalSegI](0,2) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].xz()*L()[globalSegI];

                        // blockDiag[globalSegI](1,0) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].yx()*L()[globalSegI];
                        // blockDiag[globalSegI](1,1) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].yy()*L()[globalSegI];
                        // blockDiag[globalSegI](1,2) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].yz()*L()[globalSegI];

                        // blockDiag[globalSegI](2,0) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].zx()*L()[globalSegI];
                        // blockDiag[globalSegI](2,1) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].zy()*L()[globalSegI];
                        // blockDiag[globalSegI](2,2) +=
                        //     curFrictionalContactForceDerivativeOverFcn[cI].zz()*L()[globalSegI];
                        // }

                        // //- Contact moment direction derivative
                        // if (false)
                        // // if (cc[globalSegI].size() == 2)
                        // {
                        //     const label upperNeiFace =
                        //         upperNeiFaces[globalSegI];
                        //     const label lowerNeiFace =
                        //         lowerNeiFaces[globalSegI];
                            
                        //     if (upperNeiFace < 0 || lowerNeiFace < 0)
                        //     {
                        //         FatalErrorIn
                        //         (
                        //             "coupledTotalLagNewtonRaphsonBeam::"
                        //             "applyConicalPulleysContact(...)"
                        //         )
                        //             << "Can not find upper and lower neighbour."
                        //             << abort(FatalError);
                        //     }

                        //     tensor curCoeff = 
                        //         tensor::I*(curContactMoment[cI] & T)/2;
                        //     // already multiplied by L()[globalSegI]

                        //     upper[upperNeiFace](3,0) += curCoeff.xx();
                        //     upper[upperNeiFace](3,1) += curCoeff.xy();
                        //     upper[upperNeiFace](3,2) += curCoeff.xz();

                        //     upper[upperNeiFace](4,0) += curCoeff.yx();
                        //     upper[upperNeiFace](4,1) += curCoeff.yy();
                        //     upper[upperNeiFace](4,2) += curCoeff.yz();

                        //     upper[upperNeiFace](5,0) += curCoeff.zx();
                        //     upper[upperNeiFace](5,1) += curCoeff.zy();
                        //     upper[upperNeiFace](5,2) += curCoeff.zz();

                        //     lower[lowerNeiFace](3,0) -= curCoeff.xx();
                        //     lower[lowerNeiFace](3,1) -= curCoeff.xy();
                        //     lower[lowerNeiFace](3,2) -= curCoeff.xz();

                        //     lower[lowerNeiFace](4,0) -= curCoeff.yx();
                        //     lower[lowerNeiFace](4,1) -= curCoeff.yy();
                        //     lower[lowerNeiFace](4,2) -= curCoeff.yz();

                        //     lower[lowerNeiFace](5,0) -= curCoeff.zx();
                        //     lower[lowerNeiFace](5,1) -= curCoeff.zy();
                        //     lower[lowerNeiFace](5,2) -= curCoeff.zz();
                        // }

                        // if (false)
                        // {
                        //     tensor curCoeff =
                        //        -(curContactMoment[cI] & T) * spinTensor(T);

                        //     blockDiag[globalSegI](3,3) +=
                        //         curCoeff.xx()*L()[globalSegI];
                        //     blockDiag[globalSegI](3,4) +=
                        //         curCoeff.xy()*L()[globalSegI];
                        //     blockDiag[globalSegI](3,5) +=
                        //         curCoeff.xz()*L()[globalSegI];

                        //     blockDiag[globalSegI](4,3) +=
                        //         curCoeff.yx()*L()[globalSegI];
                        //     blockDiag[globalSegI](4,4) +=
                        //         curCoeff.yy()*L()[globalSegI];
                        //     blockDiag[globalSegI](4,5) +=
                        //         curCoeff.yz()*L()[globalSegI];

                        //     blockDiag[globalSegI](5,3) +=
                        //         curCoeff.zx()*L()[globalSegI];
                        //     blockDiag[globalSegI](5,4) +=
                        //         curCoeff.zy()*L()[globalSegI];
                        //     blockDiag[globalSegI](5,5) +=
                        //         curCoeff.zz()*L()[globalSegI];
                        // }


                        
                        //- Contact moment direction derivative
                        if (false)
                        // if (cc[globalSegI].size() == 2)
                        {
                            const label upperNeiFace =
                                upperNeiFaces[globalSegI];
                            const label lowerNeiFace =
                                lowerNeiFaces[globalSegI];

                            if (upperNeiFace < 0 || lowerNeiFace < 0)
                            {
                                FatalErrorIn
                                (
                                    "coupledTotalLagNewtonRaphsonBeam::"
                                    "applyConicalPulleysContactNew(...)"
                                )
                                    << "Can not find upper and lower neighbour."
                                    << abort(FatalError);
                            }

                            tensor curCoeff = 
                                tensor::I*(curAxialContactForce[cI] & T)
                               *L()[globalSegI]*dc[upperNeiFace];
                            // already multiplied by L()[globalSegI]

                            upper[upperNeiFace](0,0) += curCoeff.xx();
                            upper[upperNeiFace](0,1) += curCoeff.xy();
                            upper[upperNeiFace](0,2) += curCoeff.xz();

                            upper[upperNeiFace](1,0) += curCoeff.yx();
                            upper[upperNeiFace](1,1) += curCoeff.yy();
                            upper[upperNeiFace](1,2) += curCoeff.yz();

                            upper[upperNeiFace](2,0) += curCoeff.zx();
                            upper[upperNeiFace](2,1) += curCoeff.zy();
                            upper[upperNeiFace](2,2) += curCoeff.zz();

                            // lower[lowerNeiFace](0,0) -= curCoeff.xx();
                            // lower[lowerNeiFace](0,1) -= curCoeff.xy();
                            // lower[lowerNeiFace](0,2) -= curCoeff.xz();

                            // lower[lowerNeiFace](1,0) -= curCoeff.yx();
                            // lower[lowerNeiFace](1,1) -= curCoeff.yy();
                            // lower[lowerNeiFace](1,2) -= curCoeff.yz();

                            // lower[lowerNeiFace](2,0) -= curCoeff.zx();
                            // lower[lowerNeiFace](2,1) -= curCoeff.zy();
                            // lower[lowerNeiFace](2,2) -= curCoeff.zz();

                            blockDiag[globalSegI](0,0) -= curCoeff.xx();
                            blockDiag[globalSegI](0,1) -= curCoeff.xy();
                            blockDiag[globalSegI](0,2) -= curCoeff.xz();

                            blockDiag[globalSegI](1,0) -= curCoeff.yx();
                            blockDiag[globalSegI](1,1) -= curCoeff.yy();
                            blockDiag[globalSegI](1,2) -= curCoeff.yz();

                            blockDiag[globalSegI](2,0) -= curCoeff.zx();
                            blockDiag[globalSegI](2,1) -= curCoeff.zy();
                            blockDiag[globalSegI](2,2) -= curCoeff.zz();
                        }
                    }
                }
            }

            start += nSeg;
        }

        // if (isA<lagrangeMultipliers>(contact().normalModel()))
        // {
        //     #include "applyLagrangeMultiplierConicalPulleyContact.H"
        // }
    }
}



void coupledTotalLagNewtonRaphsonBeam::applyToroidalPulleysContact
(
    fvBlockMatrix<vector6>& eqn
)
{
    if (toroidalPulleys().size())
    {
    	Info << "Inside toroidal pulleys contact " << endl;
        // Get matrix diagonal
        tensor6Field& blockDiag = eqn.diag().asSquare();

        // Grab source
        vector6Field& source = eqn.source();

        // // Grab off-diagonal 
        // tensor6Field& upper = eqn.upper().asSquare();

        // // Grab off-diagonal
        // tensor6Field& lower = eqn.lower().asSquare();
    
        label nPulleys = toroidalPulleys().size();

        label start = 0;
    
        label nBeams = contact().splines().size();

        for (label bI=0; bI<nBeams; bI++)
        {
            const toroidalPulleyContactListList& curPulleyContacts =
                contact().toroidalPulleyContacts()[bI];

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
            
            for
            (
                label segI=0;
                segI<contact().splines()[bI].nSegments();
                segI++
            )
            {
                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    // for (label cI=0; cI<2; cI++)
                    {
                        curq[segI] +=
                            curPulleyContacts[segI][pulleyI]
                           .normalContactForce();

                        // curq[segI] +=
                        //     curPulleyContacts[segI][pulleyI]
                        //    .frictionalContactForce();

                        // // Contact force due to friction in pulley bearings
                        // curq[segI] +=
                        //     curPulleyContacts[segI][pulleyI]
                        //    .axialFrictionalContactForce();
                    
                        // curm[segI] +=
                        //     curPulleyContacts[segI][pulleyI]
                        //    .frictionalContactMoment();
                    }
                }
            }

            forAll(curq, segI)
            {
                label globalSegI = start + segI;

                // W equation
                source[globalSegI](0) -= curq[segI].x()*L()[globalSegI];
                source[globalSegI](1) -= curq[segI].y()*L()[globalSegI];
                source[globalSegI](2) -= curq[segI].z()*L()[globalSegI];

                // Theta equation
                source[globalSegI](3) -= curm[segI].x()*L()[globalSegI];
                source[globalSegI](4) -= curm[segI].y()*L()[globalSegI];
                source[globalSegI](5) -= curm[segI].z()*L()[globalSegI];

                // b(index++) -= curm[segI].x()*L[start+segI];
                // b(index++) -= curm[segI].y()*L[start+segI];
                // b(index++) -= curm[segI].z()*L[start+segI];
                // index++;
                // index++;
                // index++;
            }

            start += curq.size();
        }

        // Info << "Add implicit contribution of conical pulley contac force"
        //      << endl;
        
        // const labelList& upperNeiFaces = upperNeiCellFaces();
        // const labelList& lowerNeiFaces = lowerNeiCellFaces();
    
        start = 0;
        for (label bI=0; bI<nBeams; bI++)
        {
            const toroidalPulleyContactListList& curPulleyContacts =
                contact().toroidalPulleyContacts()[bI];

            // const vectorField& dRdS = contact().splines()[bI].dRdS();

            label nSeg = contact().splines()[bI].nSegments();

            for (label segI=0; segI<nSeg; segI++)
            {
                label globalSegI = start + segI;

                vector T = contact().splines()[bI].dRdS()[segI];
                T /= mag(T) + SMALL;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    // const scalarField& curContactDistance =
                    //     curPulleyContacts[segI][pulleyI]
                    //    .delta();

                    // const vectorField& curContactForce =
                    //     curPulleyContacts[segI][pulleyI]
                    //    .normalContactForce();

                    const tensor& curContactForceDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForceDerivative();

                    const tensor& curContactForceDirectionDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .normalContactForceDirectionDerivative();

                    tensor FcDn =
                        curContactForceDirectionDerivative
                       *L()[globalSegI];
                    
                    // const vector& curContactMoment =
                    //     curPulleyContacts[segI][pulleyI]
                    //    .frictionalContactMoment();

                    const tensor& curContactMomentDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactMomentDerivative();

                    const tensor& curContactMomentDerivativeOverFcn =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactMomentDerivativeOverFcn();

                    const tensor& curFrictionalContactForceDerivative =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactForceDerivative();

                    const tensor&
                    curFrictionalContactForceDerivativeOverFcn =
                        curPulleyContacts[segI][pulleyI]
                       .frictionalContactForceDerivativeOverFcn();

                    // for (label cI=0; cI<2; cI++)
                    {
                        //- Contact force derivative (W)
                        blockDiag[globalSegI](0,0) +=
                            curContactForceDerivative.xx()*L()[globalSegI]
                          + FcDn.xx();
                        blockDiag[globalSegI](0,1) +=
                            curContactForceDerivative.xy()*L()[globalSegI]
                          + FcDn.xy();
                        blockDiag[globalSegI](0,2) +=
                            curContactForceDerivative.xz()*L()[globalSegI]
                          + FcDn.xz();

                        blockDiag[globalSegI](1,0) +=
                            curContactForceDerivative.yx()*L()[globalSegI]
                          + FcDn.yx();
                        blockDiag[globalSegI](1,1) +=
                            curContactForceDerivative.yy()*L()[globalSegI]
                          + FcDn.yy();
                        blockDiag[globalSegI](1,2) +=
                            curContactForceDerivative.yz()*L()[globalSegI]
                          + FcDn.yz();

                        blockDiag[globalSegI](2,0) +=
                            curContactForceDerivative.zx()*L()[globalSegI]
                          + FcDn.zx();
                        blockDiag[globalSegI](2,1) +=
                            curContactForceDerivative.zy()*L()[globalSegI]
                          + FcDn.zy();
                        blockDiag[globalSegI](2,2) +=
                            curContactForceDerivative.zz()*L()[globalSegI]
                          + FcDn.zz();
                        
                        if (false){
                        //- Frictional contact moment derivative (Theta)
                        blockDiag[globalSegI](3,3) +=
                            curContactMomentDerivative.xx()*L()[globalSegI];
                        blockDiag[globalSegI](3,4) +=
                            curContactMomentDerivative.xy()*L()[globalSegI];
                        blockDiag[globalSegI](3,5) +=
                            curContactMomentDerivative.xz()*L()[globalSegI];

                        blockDiag[globalSegI](4,3) +=
                            curContactMomentDerivative.yx()*L()[globalSegI];
                        blockDiag[globalSegI](4,4) +=
                            curContactMomentDerivative.yy()*L()[globalSegI];
                        blockDiag[globalSegI](4,5) +=
                            curContactMomentDerivative.yz()*L()[globalSegI];

                        blockDiag[globalSegI](5,3) +=
                            curContactMomentDerivative.zx()*L()[globalSegI];
                        blockDiag[globalSegI](5,4) +=
                            curContactMomentDerivative.zy()*L()[globalSegI];
                        blockDiag[globalSegI](5,5) +=
                            curContactMomentDerivative.zz()*L()[globalSegI];

                        //- Frictional contact moment derivative (W)
                        blockDiag[globalSegI](3,0) +=
                            curContactMomentDerivativeOverFcn.xx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](3,1) +=
                            curContactMomentDerivativeOverFcn.xy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](3,2) +=
                            curContactMomentDerivativeOverFcn.xz()
                           *L()[globalSegI];

                        blockDiag[globalSegI](4,0) +=
                            curContactMomentDerivativeOverFcn.yx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](4,1) +=
                            curContactMomentDerivativeOverFcn.yy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](4,2) +=
                            curContactMomentDerivativeOverFcn.yz()
                           *L()[globalSegI];

                        blockDiag[globalSegI](5,0) +=
                            curContactMomentDerivativeOverFcn.zx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](5,1) +=
                            curContactMomentDerivativeOverFcn.zy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](5,2) +=
                            curContactMomentDerivativeOverFcn.zz()
                           *L()[globalSegI];

                        //- Frictional contact force derivative (Theta)
                        if (false)
                        {
                        blockDiag[globalSegI](0,3) +=
                            curFrictionalContactForceDerivative.xx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](0,4) +=
                            curFrictionalContactForceDerivative.xy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](0,5) +=
                            curFrictionalContactForceDerivative.xz()
                           *L()[globalSegI];

                        blockDiag[globalSegI](1,3) +=
                            curFrictionalContactForceDerivative.yx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](1,4) +=
                            curFrictionalContactForceDerivative.yy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](1,5) +=
                            curFrictionalContactForceDerivative.yz()
                           *L()[globalSegI];

                        blockDiag[globalSegI](2,3) +=
                            curFrictionalContactForceDerivative.zx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](2,4) +=
                            curFrictionalContactForceDerivative.zy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](2,5) +=
                            curFrictionalContactForceDerivative.zz()
                           *L()[globalSegI];
                        }

                        //- Frictional contact force derivative (W)
                        if (false)
                        {
                        blockDiag[globalSegI](0,0) +=
                            curFrictionalContactForceDerivativeOverFcn.xx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](0,1) +=
                            curFrictionalContactForceDerivativeOverFcn.xy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](0,2) +=
                            curFrictionalContactForceDerivativeOverFcn.xz()
                           *L()[globalSegI];

                        blockDiag[globalSegI](1,0) +=
                            curFrictionalContactForceDerivativeOverFcn.yx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](1,1) +=
                            curFrictionalContactForceDerivativeOverFcn.yy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](1,2) +=
                            curFrictionalContactForceDerivativeOverFcn.yz()
                           *L()[globalSegI];

                        blockDiag[globalSegI](2,0) +=
                            curFrictionalContactForceDerivativeOverFcn.zx()
                           *L()[globalSegI];
                        blockDiag[globalSegI](2,1) +=
                            curFrictionalContactForceDerivativeOverFcn.zy()
                           *L()[globalSegI];
                        blockDiag[globalSegI](2,2) +=
                            curFrictionalContactForceDerivativeOverFcn.zz()
                           *L()[globalSegI];
                        }}
                    }
                }
            }
            
            start += nSeg;
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
