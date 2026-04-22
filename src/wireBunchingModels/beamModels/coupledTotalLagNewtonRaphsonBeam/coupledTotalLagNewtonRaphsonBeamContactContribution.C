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


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void coupledTotalLagNewtonRaphsonBeam::appendBlockTriplets
(
    std::vector<Eigen::Triplet<scalar>>& triplets,
    const label blockRow,
    const label blockCol,
    const scalarSquareMatrix& block
)
{
    for (label i = 0; i < 6; ++i)
    {
        for (label j = 0; j < 6; ++j)
        {
            triplets.emplace_back(6*blockRow + i, 6*blockCol + j, block(i,j));
        }
    }
}

void coupledTotalLagNewtonRaphsonBeam::lineContactContribution
(
    std::vector<Eigen::Triplet<scalar>>& lineContactTriplets,
    // Field<scalarSquareMatrix>& lcDiag,
    Field<scalarRectangularMatrix>& lcSource
    // PtrList<lcCoeffPairListList>& lcInterBeamCoeffs
)
{
    if (mesh().cellZones().size() > 1)
    {
        label nBeams = contact().splines().size();

        // // Get matrix diagonal
        // tensor6Field& diag = eqn.diag().asSquare();

        // // Grab source
        // vector6Field& source = eqn.source();

        // Grab beam coupling coeffs (foam-extend version)
        //        PtrList<tensor6PairListList>& lcInterBeamCoeffs = eqn.lineContactCoeffs();

        // Explicit part of the line contact force
        label start = 0;

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
                            // //-SB added
                            // const label neiSegI =
                            //                 curLineContacts[segI][nbI].secondBeamSegment();

                            // const scalar neiZeta =
                            // curLineContacts[segI][nbI].secondBeamZeta();

                            // label globalNeiSegI = whichCell(nbI, neiSegI);
                            // //-

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

                // globalSegIndex and localCellIndex are the same for serial run
                label cellID =
                    localCellIndex(globalSegIndex);

                if (cellID != -1)
                {
                    // W equation
                    lcSource[cellID](0, 0) -= curq[segI].x()*L()[cellID];
                    lcSource[cellID](1, 0) -= curq[segI].y()*L()[cellID];
                    lcSource[cellID](2, 0) -= curq[segI].z()*L()[cellID];

                    // Theta equation
                    lcSource[cellID](3, 0) -= curm[segI].x()*L()[cellID];
                    lcSource[cellID](4, 0) -= curm[segI].y()*L()[cellID];
                    lcSource[cellID](5, 0) -= curm[segI].z()*L()[cellID];
                }
            }

            start += curq.size();
        }

        // Implicit part of the line contact force
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

                                const label neiLowerSegI =
                                    curLineContacts[segI][nbI]
                                   .secondBeamLowerSegment();
                                const label neiUpperSegI =
                                    curLineContacts[segI][nbI]
                                    .secondBeamUpperSegment();

                               const label globalNeiLowerSegI =
                                   whichCell
                                   (
                                       nbI,
                                       neiLowerSegI
                                   );

                               const label globalNeiUpperSegI =
                                   whichCell
                                   (
                                       nbI,
                                       neiUpperSegI
                                   );


                                //   const scalar& curContactDistance =
                                    //   curLineContacts[segI][nbI].delta();

                                tensor curContactForceDerivative =
                                    curLineContacts[segI][nbI]
                                   .normalContactForceDerivative()
                                   *L()[cellID];

                                scalarSquareMatrix localLcDiag(6, 0.0);

                                //- Contact force derivative (W)
                                localLcDiag(0,0) =
                                    curContactForceDerivative.xx();
                                localLcDiag(0,1) =
                                    curContactForceDerivative.xy();
                                localLcDiag(0,2) =
                                    curContactForceDerivative.xz();

                                localLcDiag(1,0) =
                                    curContactForceDerivative.yx();
                                localLcDiag(1,1) =
                                    curContactForceDerivative.yy();
                                localLcDiag(1,2) =
                                    curContactForceDerivative.yz();

                                localLcDiag(2,0) =
                                    curContactForceDerivative.zx();
                                localLcDiag(2,1) =
                                    curContactForceDerivative.zy();
                                localLcDiag(2,2) =
                                    curContactForceDerivative.zz();

                                // lcDiag[cellID] += localLcDiag;


                                appendBlockTriplets(lineContactTriplets, cellID, cellID, localLcDiag);

                                // Off-diagonal

                                scalar w0 =
                                    curLineContacts[segI][nbI].secondBeamWeight();
                                scalar w1 = 1.0 - w0;

                                scalarSquareMatrix localLcInterBeamCoeffsFirst(6, 0.0);
                                scalarSquareMatrix localLcInterBeamCoeffsSecond(6, 0.0);
                                // Contact force derivative (W)

                                //- lower
                                localLcInterBeamCoeffsFirst(0, 0) =
                                    -curContactForceDerivative.xx()*w0;
                                localLcInterBeamCoeffsFirst(0, 1) =
                                    -curContactForceDerivative.xy()*w0;
                                localLcInterBeamCoeffsFirst(0, 2) =
                                    -curContactForceDerivative.xz()*w0;

                                localLcInterBeamCoeffsFirst(1, 0) =
                                    -curContactForceDerivative.yx()*w0;
                                localLcInterBeamCoeffsFirst(1, 1) =
                                    -curContactForceDerivative.yy()*w0;
                                localLcInterBeamCoeffsFirst(1, 2) =
                                    -curContactForceDerivative.yz()*w0;

                                localLcInterBeamCoeffsFirst(2, 0) =
                                    -curContactForceDerivative.zx()*w0;
                                localLcInterBeamCoeffsFirst(2, 1) =
                                    -curContactForceDerivative.zy()*w0;
                                localLcInterBeamCoeffsFirst(2, 2) =
                                    -curContactForceDerivative.zz()*w0;

                                //- upper
                                localLcInterBeamCoeffsSecond(0, 0) =
                                    -curContactForceDerivative.xx()*w1;
                                localLcInterBeamCoeffsSecond(0, 1) =
                                    -curContactForceDerivative.xy()*w1;
                                localLcInterBeamCoeffsSecond(0, 2) =
                                    -curContactForceDerivative.xz()*w1;

                                localLcInterBeamCoeffsSecond(1, 0) =
                                    -curContactForceDerivative.yx()*w1;
                                localLcInterBeamCoeffsSecond(1, 1) =
                                    -curContactForceDerivative.yy()*w1;
                                localLcInterBeamCoeffsSecond(1, 2) =
                                    -curContactForceDerivative.yz()*w1;

                                localLcInterBeamCoeffsSecond(2, 0) =
                                    -curContactForceDerivative.zx()*w1;
                                localLcInterBeamCoeffsSecond(2, 1) =
                                    -curContactForceDerivative.zy()*w1;
                                localLcInterBeamCoeffsSecond(2, 2) =
                                    -curContactForceDerivative.zz()*w1;

                                // lcInterBeamCoeffs[bI][segI][nbI].first() += localLcInterBeamCoeffsFirst;
                                // lcInterBeamCoeffs[bI][segI][nbI].second() += localLcInterBeamCoeffsSecond;

                                appendBlockTriplets
                                (
                                    lineContactTriplets,
                                    cellID, //row
                                    globalNeiLowerSegI, //column
                                    localLcInterBeamCoeffsFirst
                                );

                                appendBlockTriplets
                                (
                                    lineContactTriplets,
                                    cellID, //row
                                    globalNeiUpperSegI, //column
                                    localLcInterBeamCoeffsSecond
                                );
                            }
                        }
                    }
                }
	    }
	    start += nSeg;
	}
    }
}



void coupledTotalLagNewtonRaphsonBeam::pointContactContribution
(
    std::vector<Eigen::Triplet<scalar>>& pointContactTriplets,
    // Field<scalarSquareMatrix>& pcDiag,
    Field<scalarRectangularMatrix>& pcSource
    // Field<scalarSquareMatrix>& pcNeiUpperCoeffs,
    // Field<scalarSquareMatrix>& pcNeiLowerCoeffs,
    // pcCoeffPairListList& pcInterBeamCoeffs
)
{
    if
    (
        (mesh().cellZones().size() > 1)
     && !Pstream::parRun()
    )
    {
        // std::vector<Eigen::Triplet<scalar>> pointContactTriplets;
        // pointContactTriplets.reserve(300*contact().pointContacts().size());
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
    //    tensor6PairListList& pcInterBeamCoeffs = eqn.pointContactCoeffs();


    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();

    const surfaceVectorField dRdS(dR0Ds_ + fvc::snGrad(W_));

    const surfaceScalarField& dc = mesh().deltaCoeffs();


    forAll(contact().pointContacts(), pcI)
    {

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
            globalSegI[0] = whichCell(bI[0], segI);

	    // SB (2023 May): changed the explicit contact contribution
            contactForce[0] =
                contact().pointContacts()[pcI].normalContactForce();
    	    tangContactForce[0] =
                contact().pointContacts()[pcI].firstBeamTangContactForce()
              - contact().pointContacts()[pcI].secondBeamTangContactForce();

        }

        // Neighbour
        {
            bI[1] = contact().pointContacts()[pcI].secondBeam();
            label neiSegI = contact().pointContacts()[pcI].secondBeamSegment();
            zeta[1] = contact().pointContacts()[pcI].secondBeamZeta();
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

            Info<< "Tangent contact Force " << tangContactForce[sI] << endl;
            Info<< "Mag normal contact force " << mag(contactForce[sI]) << endl;

            vector F0 = contactForce[sI] + tangContactForce[sI];

            vector M0 =  (spinTensor(DR) & F0);

            // W equation
            pcSource[globalSegI[sI]](0, 0) -= F0.x();
            pcSource[globalSegI[sI]](1, 0) -= F0.y();
            pcSource[globalSegI[sI]](2, 0) -= F0.z();


            // Theta equation
            pcSource[globalSegI[sI]](3, 0) -= M0.x();
            pcSource[globalSegI[sI]](4, 0) -= M0.y();
            pcSource[globalSegI[sI]](5, 0) -= M0.z();

            Info<< "pcSource: " <<  pcSource[globalSegI[sI]](1, 0) << endl;
        }
    }


    // forAll(pcSource, cellI)
    //     {
    //         Info<< "cellI: " << cellI <<  "pcSource: " << pcSource[cellI] << endl;
    //     };
    // // Create Eigen sparse matrix and set coeffs
    // Eigen::SparseMatrix<scalar> A;

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

        // Info<<"global Segments: " << globalSegI[0] << tab << globalSegI[1]
        //     << "\nglobal lower seg: " << globalLowerSegI[0] << tab << globalLowerSegI[1]
        //     << "\nglobal Upper Segment: " << globalUpperSegI[0] << tab << globalUpperSegI[1]
        //     << endl;

        forAll(globalSegI, sI)
        {
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

            Info<< "wDiag for global seg: " << globalSegI[sI] << " is " << wDiag << endl;
            Info<< "For beam " << bI[0] << " - Shared face between "
                << "cells: " << globalLowerSegI[sI] << " and "
                << globalUpperSegI[sI] << " is face "
                << sharedFaceIndex << " and owner of face "
                << sharedFaceIndex << " is " << own[sharedFaceIndex]
                << endl;

            // Info<< "diagCell: " << diagCell << " , and wDiag " << wDiag << endl;
            Info<< "contactForceDerivative:\n" << contactForceDerivative[sI]
                << "\ncontactForceDn:\n" << contactForceDn[sI]
                << endl;

            scalarSquareMatrix localPcDiag(6, 0.0);

            //-- wDiag
            //- Force
            localPcDiag(0, 0) =
                wDiag*contactForceDerivative[sI].xx()
              + wDiag*tangContactForceDerivative[sI].xx()
    	      + wDiag*contactForceDn[sI].xx()
	          + wDiag*tangContactForceDt[sI].xx();

            localPcDiag(0, 1) =
                wDiag*contactForceDerivative[sI].xy()
              + wDiag*tangContactForceDerivative[sI].xy()
    	      + wDiag*contactForceDn[sI].xy()
              + wDiag*tangContactForceDt[sI].xy();

            localPcDiag(0, 2) =
                wDiag*contactForceDerivative[sI].xz()
              + wDiag*tangContactForceDerivative[sI].xz()
    	      + wDiag*contactForceDn[sI].xz()
              + wDiag*tangContactForceDt[sI].xz();

            localPcDiag(1, 0) =
                wDiag*contactForceDerivative[sI].yx()
              + wDiag*tangContactForceDerivative[sI].yx()
    	      + wDiag*contactForceDn[sI].yx()
	          + wDiag*tangContactForceDt[sI].yx();

            localPcDiag(1, 1) =
                wDiag*contactForceDerivative[sI].yy()
              + wDiag*tangContactForceDerivative[sI].yy()
    	      + wDiag*contactForceDn[sI].yy()
	          + wDiag*tangContactForceDt[sI].yy();

            localPcDiag(1, 2) =
                wDiag*contactForceDerivative[sI].yz()
              + wDiag*tangContactForceDerivative[sI].yz()
    	      + wDiag*contactForceDn[sI].yz()
	          + wDiag*tangContactForceDt[sI].yz();


            localPcDiag(2, 0) =
                wDiag*contactForceDerivative[sI].zx()
              + wDiag*tangContactForceDerivative[sI].zx()
    	      + wDiag*contactForceDn[sI].zx()
	          + wDiag*tangContactForceDt[sI].zx();

            localPcDiag(2, 1) =
                wDiag*contactForceDerivative[sI].zy()
              + wDiag*tangContactForceDerivative[sI].zy()
              + wDiag*contactForceDn[sI].zy()
    	      + wDiag*tangContactForceDt[sI].zy();

            localPcDiag(2, 2) =
                wDiag*contactForceDerivative[sI].zz()
              + wDiag*tangContactForceDerivative[sI].zz()
              + wDiag*contactForceDn[sI].zz()
    	      + wDiag*tangContactForceDt[sI].zz();


            //- Moment
            localPcDiag(3, 0) = wDiag*DMCoeff.xx();
            localPcDiag(3, 1) = wDiag*DMCoeff.xy();
            localPcDiag(3, 2) = wDiag*DMCoeff.xz();

            localPcDiag(4, 0) = wDiag*DMCoeff.yx();
            localPcDiag(4, 1) = wDiag*DMCoeff.yy();
            localPcDiag(4, 2) = wDiag*DMCoeff.yz();

            localPcDiag(5, 0) = wDiag*DMCoeff.zx();
            localPcDiag(5, 1) = wDiag*DMCoeff.zy();
            localPcDiag(5, 2) = wDiag*DMCoeff.zz();


            appendBlockTriplets(pointContactTriplets, diagCell, diagCell, localPcDiag);

            // pcDiag[diagCell] += localPcDiag;

            //-- w1
            if (sharedFaceIndex != -1)
            {
                scalarSquareMatrix curCoeff(6, 0.0);

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


                if (own[sharedFaceIndex] != globalSegI[sI])
                {
                    // pcNeiLowerCoeffs[sharedFaceIndex] += curCoeff;

                    appendBlockTriplets
                    (
                        pointContactTriplets,
                        nei[sharedFaceIndex], // row
                        own[sharedFaceIndex], // column
                        curCoeff
                    );
                }
                else
                {
                    // pcNeiUpperCoeffs[sharedFaceIndex] += curCoeff;

                    appendBlockTriplets
                    (
                        pointContactTriplets,
                        own[sharedFaceIndex], // row
                        nei[sharedFaceIndex], // column
                        curCoeff
                    );
                }
            }


            // Off-diagonal
            // label globalNeiSeg0 = globalLowerNeiSegI[sI];
            // label globalNeiSeg1 = globalUpperNeiSegI[sI];
            // label globalNeiSeg0 = globalNeiSegI[sI];
            // label globalNeiSeg1 = globalNeiSegI[sI];

            scalar w0 = neiWeight[sI];
            scalar w1 = 1.0 - w0;

            Info<< "w0: " << w0 << endl;

            scalarSquareMatrix pcLocalInterBeamCoeffsFirst(6, 0.0);
            scalarSquareMatrix pcLocalInterBeamCoeffsSecond(6, 0.0);

            //-- IN0
            //- Force
            pcLocalInterBeamCoeffsFirst(0,0) -=
                w0*contactForceDerivative[sI].xx()
              + w0*contactForceDn[sI].xx()
    	      + w0*couplingTangentContactForceCoeffs[sI].xx();
            pcLocalInterBeamCoeffsFirst(0,1) -=
                w0*contactForceDerivative[sI].xy()
              + w0*contactForceDn[sI].xy()
	          + w0*couplingTangentContactForceCoeffs[sI].xy();
            pcLocalInterBeamCoeffsFirst(0,2) -=
                w0*contactForceDerivative[sI].xz()
              + w0*contactForceDn[sI].xz()
	          + w0*couplingTangentContactForceCoeffs[sI].xz();

            pcLocalInterBeamCoeffsFirst(1,0) -=
                w0*contactForceDerivative[sI].yx()
              + w0*contactForceDn[sI].yx()
	          + w0*couplingTangentContactForceCoeffs[sI].yx();
            pcLocalInterBeamCoeffsFirst(1,1) -=
                w0*contactForceDerivative[sI].yy()
              + w0*contactForceDn[sI].yy()
    	      + w0*couplingTangentContactForceCoeffs[sI].yy();
            pcLocalInterBeamCoeffsFirst(1,2) -=
                w0*contactForceDerivative[sI].yz()
              + w0*contactForceDn[sI].yz()
	          + w0*couplingTangentContactForceCoeffs[sI].yz();

            pcLocalInterBeamCoeffsFirst(2,0) -=
                w0*contactForceDerivative[sI].zx()
              + w0*contactForceDn[sI].zx()
	          + w0*couplingTangentContactForceCoeffs[sI].zx();
            pcLocalInterBeamCoeffsFirst(2,1) -=
                w0*contactForceDerivative[sI].zy()
              + w0*contactForceDn[sI].zy()
    	      + w0*couplingTangentContactForceCoeffs[sI].zy();
            pcLocalInterBeamCoeffsFirst(2,2) -=
                w0*contactForceDerivative[sI].zz()
              + w0*contactForceDn[sI].zz()
    	      + w0*couplingTangentContactForceCoeffs[sI].zz();


            //- Moment
            pcLocalInterBeamCoeffsFirst(3,0) -=
                w0*DMCoeff.xx();
            pcLocalInterBeamCoeffsFirst(3,1) -=
                w0*DMCoeff.xy();
            pcLocalInterBeamCoeffsFirst(3,2) -=
                w0*DMCoeff.xz();

            pcLocalInterBeamCoeffsFirst(4,0) -=
                w0*DMCoeff.yx();
            pcLocalInterBeamCoeffsFirst(4,1) -=
                w0*DMCoeff.yy();
            pcLocalInterBeamCoeffsFirst(4,2) -=
                w0*DMCoeff.yz();

            pcLocalInterBeamCoeffsFirst(5,0) -=
                w0*DMCoeff.zx();
            pcLocalInterBeamCoeffsFirst(5,1) -=
                w0*DMCoeff.zy();
            pcLocalInterBeamCoeffsFirst(5,2) -=
                w0*DMCoeff.zz();

            //-- IN1
            //- Force
            pcLocalInterBeamCoeffsSecond(0,0) -=
                w1*contactForceDerivative[sI].xx()
              + w1*contactForceDn[sI].xx()
	          + w1*couplingTangentContactForceCoeffs[sI].xx();
            pcLocalInterBeamCoeffsSecond(0,1) -=
                w1*contactForceDerivative[sI].xy()
              + w1*contactForceDn[sI].xy()
	          + w1*couplingTangentContactForceCoeffs[sI].xy();
            pcLocalInterBeamCoeffsSecond(0,2) -=
                w1*contactForceDerivative[sI].xz()
              + w1*contactForceDn[sI].xz()
    	      + w1*couplingTangentContactForceCoeffs[sI].xz();

            pcLocalInterBeamCoeffsSecond(1,0) -=
                w1*contactForceDerivative[sI].yx()
              + w1*contactForceDn[sI].yx()
    	      + w1*couplingTangentContactForceCoeffs[sI].yx();
            pcLocalInterBeamCoeffsSecond(1,1) -=
                w1*contactForceDerivative[sI].yy()
              + w1*contactForceDn[sI].yy()
	          + w1*couplingTangentContactForceCoeffs[sI].yy();
            pcLocalInterBeamCoeffsSecond(1,2) -=
                w1*contactForceDerivative[sI].yz()
              + w1*contactForceDn[sI].yz()
	          + w1*couplingTangentContactForceCoeffs[sI].yz();

            pcLocalInterBeamCoeffsSecond(2,0) -=
                w1*contactForceDerivative[sI].zx()
              + w1*contactForceDn[sI].zx()
    	      + w1*couplingTangentContactForceCoeffs[sI].zx();
            pcLocalInterBeamCoeffsSecond(2,1) -=
                w1*contactForceDerivative[sI].zy()
              + w1*contactForceDn[sI].zy()
	          + w1*couplingTangentContactForceCoeffs[sI].zy();
            pcLocalInterBeamCoeffsSecond(2,2) -=
                w1*contactForceDerivative[sI].zz()
              + w1*contactForceDn[sI].zz()
	          + w1*couplingTangentContactForceCoeffs[sI].zz();

            //- Moment
            pcLocalInterBeamCoeffsSecond(3,0) -=
                w1*DMCoeff.xx();
            pcLocalInterBeamCoeffsSecond(3,1) -=
                w1*DMCoeff.xy();
            pcLocalInterBeamCoeffsSecond(3,2) -=
                w1*DMCoeff.xz();

            pcLocalInterBeamCoeffsSecond(4,0) -=
                w1*DMCoeff.yx();
            pcLocalInterBeamCoeffsSecond(4,1) -=
                w1*DMCoeff.yy();
            pcLocalInterBeamCoeffsSecond(4,2) -=
                w1*DMCoeff.yz();

            pcLocalInterBeamCoeffsSecond(5,0) -=
                w1*DMCoeff.zx();
            pcLocalInterBeamCoeffsSecond(5,1) -=
                w1*DMCoeff.zy();
            pcLocalInterBeamCoeffsSecond(5,2) -=
                w1*DMCoeff.zz();

            appendBlockTriplets
            (
                pointContactTriplets,
                globalSegI[sI], // row
                globalLowerNeiSegI[sI], // column
                pcLocalInterBeamCoeffsFirst
            );

            appendBlockTriplets
            (
                pointContactTriplets,
                globalSegI[sI], // row
                globalUpperNeiSegI[sI], // column
                pcLocalInterBeamCoeffsSecond
            );

            // Info<< "Global seg I: " << globalSegI[sI]
            //     << "\nGlobal Lower NeiSeg I: " << globalLowerNeiSegI[sI]
            //     << "\nGlobal Upper NeiSeg I: " << globalUpperNeiSegI[sI]
            //     << endl;


            // pcInterBeamCoeffs[pcI][sI].first() += pcLocalInterBeamCoeffsFirst;
            // pcInterBeamCoeffs[pcI][sI].second() += pcLocalInterBeamCoeffsSecond;

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
