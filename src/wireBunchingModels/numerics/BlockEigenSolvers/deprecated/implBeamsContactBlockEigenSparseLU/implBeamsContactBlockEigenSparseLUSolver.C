/*---------------------------------------------------------------------------* \
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

\*---------------------------------------------------------------------------*/

#include "blockLduSolvers.H"
#include "implBeamsContactBlockEigenSparseLUSolver.H"
#include "addToRunTimeSelectionTable.H"
#include <Eigen/Sparse>
// #include <GMRES.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
// #include <unsupported/Eigen/IterativeSolvers>

// #include <Eigen/PardisoSupport>


#include "fvMesh.H"
#include "fvc.H"
#include "beamModel.H"
#include "cubicSpline.H"
#include "Tuple2.H"
#include "vectorList.H"
#include "spinTensor.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(implBeamsContactBlockEigenSparseLUSolver, 0);
    addToRunTimeSelectionTable
    (
        blockVector6Solver, implBeamsContactBlockEigenSparseLUSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockVector6Solver, implBeamsContactBlockEigenSparseLUSolver, asymMatrix
    );

    typedef Tuple2<label, scalar> labelScalar;
    typedef List<labelScalar> labelScalarList;
    typedef List<labelScalarList> labelScalarListList;

    typedef List<tensor> tensorList;
    typedef List<tensorList> tensorListList;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::implBeamsContactBlockEigenSparseLUSolver::
implBeamsContactBlockEigenSparseLUSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector6>& matrix,
    const dictionary& dict
)
:
    BlockEigenSolver(fieldName, matrix, dict)
{
    // Info<< type() << " Eigen linear solver" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::BlockSolverPerformance<Foam::vector6>
Foam::implBeamsContactBlockEigenSparseLUSolver::solve
(
    Field<Foam::vector6>& U,
    const Field<Foam::vector6>& blockB
)
{
    if (Pstream::parRun())
    {
        FatalErrorIn
        (
            "Foam::BlockSolverPerformance<Foam::vector>"
            "Foam::implBeamContactBlockEigenSparseLUSolver::solve"
            "("
            "    Field<Foam::vector>& U,"
            "    const Field<Foam::vector>& blockB"
            ")"
        )   << "SparesLU direct linear solver may not be run in parallel"
            << abort(FatalError);
    }

    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();

    // Calculate te number of degrees of freedom
    const int m = calcDegreesOfFreedom(matrix_, twoD);

    // Create Eigen sparse matrix and set coeffs
    Eigen::SparseMatrix<scalar> A(m, m);

    // Convert foam matrix to the Eigen matrix format
    convertFoamMatrixToEigenMatrix(matrix_, A);

    // Copy source vector into Eigen vector
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(m);
    label index = 0;
    forAll(blockB, rowI)
    {
        b(index++) = blockB[rowI](0);
        b(index++) = blockB[rowI](1);
        b(index++) = blockB[rowI](2);
        b(index++) = blockB[rowI](3);
        b(index++) = blockB[rowI](4);
        b(index++) = blockB[rowI](5);
    }

    // Add implicit contact
    addImplicitContact(A, b);

    // Copy source vector into Eigen vector
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(m);
    index = 0;
    forAll(U, rowI)
    {
        x(index++) = U[rowI](0);
        x(index++) = U[rowI](1);
        x(index++) = U[rowI](2);
        x(index++) = U[rowI](3);
        x(index++) = U[rowI](4);
        x(index++) = U[rowI](5);
    }

    // Eigen::Matrix<scalar, Eigen::Dynamic, 1> r = A*x-b;
    double ir = (A*x - b).norm() / b.norm();

    // Optionally export system to Matlab
    if (writeMatlabFiles())
    {
        writeLinearSystemToMatlabFiles(A, b);
    }

    // Compressing matrix is meant to help performance
    A.makeCompressed();

    // Create scaling
    //Eigen::IterScaling<Eigen::SparseMatrix<double> >* scalPtr = NULL;

    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": direct solution of sparse system" << endl;
    }

    // Now, solve the equilibrated linear system with any available solver
    // In this case we use a Direct sparse LU solver
    // Trial an error suggests COLAMD ordering is by far the best
    // Eigen::SimplicialLDLT
    Eigen::SparseLU
    <
        Eigen::SparseMatrix<scalar>, Eigen::COLAMDOrdering<int>
        // Eigen::SparseMatrix<scalar>, Eigen::AMDOrdering<int>
        // Eigen::SparseMatrix<scalar>, Eigen::MetisOrdering<int>
    > solver(A);

    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<scalar> > solver(A);
    // Eigen::SparseQR<Eigen::SparseMatrix<scalar>, Eigen::COLAMDOrdering<int> > solver(A);
    // x = A.colPivHouseholderQr().solve(b);

    x = solver.solve(b);
    // Eigen::Matrix<scalar, Eigen::Dynamic, 1> x = solver.solve(b);

    // if (false)
    // {
    //     Eigen::PardisoLU
    //     <
    //         Eigen::SparseMatrix<scalar>
    //     > solver(A);
        
    //     x = solver.solve(b);
    // }
    
    // {
    //     Eigen::GMRES<Eigen::SparseMatrix<scalar>, Eigen::IdentityPreconditioner> gmres;
    //     gmres.compute(A);
    //     x = gmres.solve(b);
    //     std::cout << "GMRES:    #iterations: " << gmres.iterations() << ", estimated error: " << gmres.error() << std::endl;
    // }

    
    double fr = (A*x - b).norm() / b.norm();
    
    // We copy the results from the std::vector into the geometric field
    // This can be costly
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName << ": copying results into foam format" << nl
            << endl;
    }
    index = 0;
    forAll(U, cellI)
    {
        U[cellI](0) = x(index++);
        U[cellI](1) = x(index++);
        U[cellI](2) = x(index++);
        U[cellI](3) = x(index++);
        U[cellI](4) = x(index++);
        U[cellI](5) = x(index++);

        // if (!twoD)
        // {

        //     U[cellI].z() = x(index++);
        // }
    }

    Foam::vector6 initRes = Foam::vector6::zero;
    initRes(0) = ir;
    // initRes(1) = ir;
    // initRes(2) = ir;
    // initRes(3) = ir;
    // initRes(4) = ir;
    // initRes(5) = ir;
    
    Foam::vector6 finalRes = Foam::vector6::zero;
    finalRes(0) = fr;
    // finalRes(1) = fr;
    // finalRes(2) = fr;
    // finalRes(3) = fr;
    // finalRes(4) = fr;
    // finalRes(5) = fr;
    
    return BlockSolverPerformance<Foam::vector6>
    (
        this->typeName,
        this->fieldName(),
        initRes,
        // pTraits<Foam::vector6>::one,
        finalRes,
        // pTraits<Foam::vector6>::zero,
        0,
        true,
        false
    );
}

void Foam::implBeamsContactBlockEigenSparseLUSolver::addImplicitContact
(
    Eigen::SparseMatrix<scalar>& A,
    Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b
)
{
    const fvMesh& mesh
    (
        refCast<const fvMesh>
        (
            matrix_.mesh()
        )
    );

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>
        (
            "beamProperties"
        );

    if (!beam.contactActive())
    {
        return;
    }
    
    label nBeams = beam.contact().splines().size();

    // Set explicit part of the total contact force
    const scalarField& L = beam.L();
    label start = 0;
    label index = 0;
    if (false){
    for (label bI=0; bI<nBeams; bI++)
    {
        const vectorListList& curContactForces =
            beam.contact().contactForces()[bI];

        vectorField curq
        (
            beam.contact().splines()[bI].nSegments(),
            vector::zero
        );

        forAll(curContactForces, nbI)
        {
            if (nbI!=bI) // No self contact
            {
                curq += curContactForces[nbI];
            }
        }

        forAll(curq, segI)
        {
            b(index++) -= curq[segI].x()*L[start+segI];
            b(index++) -= curq[segI].y()*L[start+segI];
            b(index++) -= curq[segI].z()*L[start+segI];
            // Skip theta equation
            index++;
            index++;
            index++;
        }

        start += curq.size();
    }
    }
    else
    {
    for (label bI=0; bI<nBeams; bI++)
    {
        const lineContactListList& curLineContacts =
	    beam.contact().lineContacts()[bI];

        vectorField curq
        (
            beam.contact().splines()[bI].nSegments(),
            vector::zero
        );
        vectorField curm
        (
            beam.contact().splines()[bI].nSegments(),
            vector::zero
        );

	for
	(
	    label segI=0;
	    segI<beam.contact().splines()[bI].nSegments();
	    segI++
	)
	{
            for (label nbI=0; nbI<nBeams; nbI++)
            {
                if (nbI != bI) // No self-contact
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

        forAll(curq, segI)
        {
            // W equation
            b(index++) -= curq[segI].x()*L[start+segI];
            b(index++) -= curq[segI].y()*L[start+segI];
            b(index++) -= curq[segI].z()*L[start+segI];
            // Theta equation
            b(index++) -= curm[segI].x()*L[start+segI];
            b(index++) -= curm[segI].y()*L[start+segI];
            b(index++) -= curm[segI].z()*L[start+segI];
            // index++;
            // index++;
            // index++;

	    // Info << bI << ", " << segI << ", " << curq[segI] << ", " << curm[segI] << endl;
        }

        start += curq.size();
    }
    }

    // Add explicit part of point contact forces
    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();
    const surfaceVectorField& dR0Ds =
        mesh.lookupObject<surfaceVectorField>("dR0Ds");
    const volVectorField& W = mesh.lookupObject<volVectorField>("W");
    const surfaceVectorField dRdS = dR0Ds + fvc::snGrad(W);    
    const surfaceScalarField& dc = mesh.deltaCoeffs();
    start = 0;
    forAll(beam.contact().pointContacts(), pcI)
    {
        labelList bI(2, -1);
        labelList globalSegI(2, -1);
        scalarField zeta(2, -2);
        vectorList contactForce(2, vector::zero);
        vectorList tangContactForce(2, vector::zero);

        //Owner
        {
            bI[0] = beam.contact().pointContacts()[pcI].firstBeam();
            label segI = beam.contact().pointContacts()[pcI].firstBeamSegment();
            zeta[0] = beam.contact().pointContacts()[pcI].firstBeamZeta();
            label start = 0;
            label i = 0;
            while(i < bI[0])
            {
                start += beam.contact().splines()[i].nSegments();
                i++;
            }
            globalSegI[0] = start + segI;
            contactForce[0] =
                beam.contact().pointContacts()[pcI].normalContactForce()
              + beam.contact().pointContacts()[pcI].firstBeamTangContactForce()
              - beam.contact().pointContacts()[pcI].secondBeamTangContactForce();

            // tangContactForce[0] =
            //   beam.contact().pointContacts()[pcI].firstBeamTangContactForce();
            
            // Info << beam.contact().pointContacts()[pcI].firstBeamTangContactForce() << endl;
            // Info << beam.contact().pointContacts()[pcI].secondBeamTangContactForce() << endl;
        }

        // Neighbour
        {
            bI[1] = beam.contact().pointContacts()[pcI].secondBeam();
            label neiSegI = beam.contact().pointContacts()[pcI].secondBeamSegment();
            zeta[1] = beam.contact().pointContacts()[pcI].secondBeamZeta();
            label neiStart = 0;
            label i = 0;
            while(i < bI[1])
            {
                neiStart += beam.contact().splines()[i].nSegments();
                i++;
            }
            globalSegI[1] = neiStart + neiSegI;
            
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
                        mesh.boundary()[beam.endPatchIndex(bI[sI])]
                       .faceCells();

                    label bFaceID = findIndex(faceCells, globalSegI[sI]);
                
                    DR = zeta[sI]*dRdS.boundaryField()[beam.endPatchIndex(bI[sI])][bFaceID]
                       /dc.boundaryField()[beam.endPatchIndex(bI[sI])][bFaceID];
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
                        mesh.boundary()[beam.startPatchIndex(bI[sI])]
                       .faceCells();

                    label bFaceID = findIndex(faceCells, globalSegI[sI]);

                    DR = zeta[sI]*dRdS.boundaryField()[beam.startPatchIndex(bI[sI])][bFaceID]
                       /dc.boundaryField()[beam.startPatchIndex(bI[sI])][bFaceID];
                }
                else
                {
                    DR = 0.5*zeta[sI]*dRdS.internalField()[faceID]
                       /dc.internalField()[faceID];
                }
            }

            vector F0 = contactForce[sI];
            vector M0 = (spinTensor(DR) & F0);
            
            label index = 6*globalSegI[sI];
            b(index++) -= F0.x();
            b(index++) -= F0.y();
            b(index++) -= F0.z();
            b(index++) -= M0.x();
            b(index++) -= M0.y();
            b(index++) -= M0.z();
        }
    }
      
    
    // Set implicit part of the total contact force
    if (beam.contact().implicit())
    {
        if (false) {
        start = 0;
        for (label bI=0; bI<nBeams; bI++)
        {
            const labelList& curBeamCells = mesh.cellZones()[bI];

            label neiStart = 0;
            for (label nbI=0; nbI<nBeams; nbI++)
            {
                if (nbI!=bI) // No self contact
                {                    
                    const scalarList& curContactDistances =
                        beam.contact().contactDistances()[bI][nbI];

                    const vectorList& curContactForces =
                        beam.contact().contactForces()[bI][nbI];

                    const tensorList& curContactForceDerivatives =
                        beam.contact().contactForceDerivatives()[bI][nbI];

                    const labelScalarList& curNearestPoints =
                        beam.contact().nearestPoints()[bI][nbI];

                    label nSegments = beam.contact().splines()[nbI].nSegments();
                    scalarField neiL =
                        beam.contact().splines()[nbI].segLengths();

                    vectorField dRdZeta =
                        beam.contact().splines()[nbI]
                       .firstDerivativeParam(curNearestPoints);
                    vectorField d2RdZeta2 =
                        beam.contact().splines()[nbI]
                       .secondDerivativeParam(curNearestPoints);
                    
                    tensorField TT(curContactForces.size(), tensor::zero);
                    tensorField ImNN(curContactForces.size(), tensor::zero);
                    forAll(ImNN, segI)
                    {
                        label neiSegI = curNearestPoints[segI].first();
                        
                        vector N = curContactForces[segI];
                        if (mag(N) > SMALL)
                        {
                            N /= mag(N);
                            ImNN[segI] = I - (N*N);
                        }
                        
                        TT[segI] = (dRdZeta[neiSegI]*dRdZeta[neiSegI])
                           /((dRdZeta[neiSegI]&dRdZeta[neiSegI]));

                        // TT[segI] = tensor::zero;
                        
                        // vector d = (curContactDistances[segI]*N);
                        // TT[segI] = (dRdZeta[neiSegI]*dRdZeta[neiSegI])
                        //    /(
                        //         (dRdZeta[neiSegI]&dRdZeta[neiSegI])
                        //       - (d & d2RdZeta2[neiSegI])
                        //     );
                    }

                    forAll(curBeamCells, segI)
                    {
                        // Diagonal
                        label globalSegI = start + segI;
                        label ID = 6*globalSegI;

                        label neiSegI = curNearestPoints[segI].first();

                        tensor FcDn =
                        (
                            mag(curContactForces[segI])
                           *(ImNN[segI] - TT[segI])
                           *L[globalSegI]/curContactDistances[segI]
                        );

                        // FcDn = tensor::zero;

                        A.coeffRef(ID,ID) +=
                            curContactForceDerivatives[segI].xx()*L[globalSegI]
                          + FcDn.xx();
                        A.coeffRef(ID,ID+1) +=
                            curContactForceDerivatives[segI].xy()*L[globalSegI]
                          + FcDn.xy();
                        A.coeffRef(ID,ID+2) +=
                            curContactForceDerivatives[segI].xz()*L[globalSegI]
                          + FcDn.xz();

                        A.coeffRef(ID+1,ID) +=
                            curContactForceDerivatives[segI].yx()*L[globalSegI]
                          + FcDn.yx();
                        A.coeffRef(ID+1,ID+1) +=
                            curContactForceDerivatives[segI].yy()*L[globalSegI]
                          + FcDn.yy();
                        A.coeffRef(ID+1,ID+2) +=
                            curContactForceDerivatives[segI].yz()*L[globalSegI]
                          + FcDn.yz();

                        A.coeffRef(ID+2,ID) +=
                            curContactForceDerivatives[segI].zx()*L[globalSegI]
                          + FcDn.zx();
                        A.coeffRef(ID+2,ID+1) +=
                            curContactForceDerivatives[segI].zy()*L[globalSegI]
                          + FcDn.zy();
                        A.coeffRef(ID+2,ID+2) +=
                            curContactForceDerivatives[segI].zz()*L[globalSegI]
                          + FcDn.zz();

                        // Off-diagonal
                        // label neiSegI = curNearestPoints[segI].first();
                        scalar neiSegParam = curNearestPoints[segI].second();
                        label globalNeiSegI = neiStart + neiSegI;
                        label globalNeiSeg0 = globalNeiSegI;
                        label globalNeiSeg1 = globalNeiSegI;
                        scalar w0 = 0.5;
                        scalar w1 = 0.5;

                        if
                        (
                            (neiSegI > 0)
                         && (neiSegI < (nSegments-1))
                        )
                        {
                            if (neiSegParam >= 0)
                            {
                                globalNeiSeg1 = globalNeiSegI + 1;
                                scalar l = 0.5*neiL[neiSegI]
                                  + 0.5*neiL[neiSegI+1];
                                scalar l0 = neiSegParam*neiL[neiSegI]/2;
                                w0 = (l-l0)/l;
                                w1 = 1-w0;
                            }
                            else
                            {
                                globalNeiSeg0 = globalNeiSegI - 1;
                                scalar l = 0.5*neiL[neiSegI-1]
                                  + 0.5*neiL[neiSegI];
                                scalar l1 = -neiSegParam*neiL[neiSegI]/2;
                                w1 = (l-l1)/l;
                                w0 = 1-w1;
                            }
                        }

                        label IN0 = 6*globalNeiSeg0;
                        label IN1 = 6*globalNeiSeg1;

                        // IN0
                        A.coeffRef(ID,IN0) -=
                            curContactForceDerivatives[segI].xx()
                           *L[globalSegI]*w0 + FcDn.xx()*w0;
                        A.coeffRef(ID,IN0+1) -=
                            curContactForceDerivatives[segI].xy()
                           *L[globalSegI]*w0 + FcDn.xy()*w0;
                        A.coeffRef(ID,IN0+2) -=
                            curContactForceDerivatives[segI].xz()
                           *L[globalSegI]*w0 + FcDn.xz()*w0;

                        A.coeffRef(ID+1,IN0) -=
                            curContactForceDerivatives[segI].yx()
                           *L[globalSegI]*w0 + FcDn.yx()*w0;
                        A.coeffRef(ID+1,IN0+1) -=
                            curContactForceDerivatives[segI].yy()
                           *L[globalSegI]*w0 + FcDn.yy()*w0;
                        A.coeffRef(ID+1,IN0+2) -=
                            curContactForceDerivatives[segI].yz()
                           *L[globalSegI]*w0 + FcDn.yz()*w0;

                        A.coeffRef(ID+2,IN0) -=
                            curContactForceDerivatives[segI].zx()
                           *L[globalSegI]*w0 + FcDn.zx()*w0;
                        A.coeffRef(ID+2,IN0+1) -=
                            curContactForceDerivatives[segI].zy()
                           *L[globalSegI]*w0 + FcDn.zy()*w0;
                        A.coeffRef(ID+2,IN0+2) -=
                            curContactForceDerivatives[segI].zz()
                           *L[globalSegI]*w0 + FcDn.zz()*w0;

                        // IN1
                        A.coeffRef(ID,IN1) -=
                            curContactForceDerivatives[segI].xx()
                           *L[globalSegI]*w1 + FcDn.xx()*w1;
                        A.coeffRef(ID,IN1+1) -=
                            curContactForceDerivatives[segI].xy()
                           *L[globalSegI]*w1 + FcDn.xy()*w1;
                        A.coeffRef(ID,IN1+2) -=
                            curContactForceDerivatives[segI].xz()
                           *L[globalSegI]*w1 + FcDn.xz()*w1;

                        A.coeffRef(ID+1,IN1) -=
                            curContactForceDerivatives[segI].yx()
                           *L[globalSegI]*w1 + FcDn.yx()*w1;
                        A.coeffRef(ID+1,IN1+1) -=
                            curContactForceDerivatives[segI].yy()
                           *L[globalSegI]*w1 + FcDn.yy()*w1;
                        A.coeffRef(ID+1,IN1+2) -=
                            curContactForceDerivatives[segI].yz()
                           *L[globalSegI]*w1 + FcDn.yz()*w1;

                        A.coeffRef(ID+2,IN1) -=
                            curContactForceDerivatives[segI].zx()
                           *L[globalSegI]*w1 + FcDn.zx()*w1;
                        A.coeffRef(ID+2,IN1+1) -=
                            curContactForceDerivatives[segI].zy()
                           *L[globalSegI]*w1 + FcDn.zy()*w1;
                        A.coeffRef(ID+2,IN1+2) -=
                            curContactForceDerivatives[segI].zz()
                           *L[globalSegI]*w1 + FcDn.zz()*w1;
                    }
                }

                neiStart += beam.contact().splines()[nbI].nSegments();
            }

            start += curBeamCells.size();
        }
	}
	else
	{
	label start = 0;
        for (label bI=0; bI<nBeams; bI++)
	{
            const lineContactListList& curLineContacts =
	        beam.contact().lineContacts()[bI];

	    label nSeg = beam.contact().splines()[bI].nSegments();

	    for (label segI=0; segI<nSeg; segI++)
	    {
                for (label nbI=0; nbI<nBeams; nbI++)
		{
                    if (nbI != bI) // No self-contact
		    {
			const label neiSegI =
  			    curLineContacts[segI][nbI].secondBeamSegment();

                        const vector& curContactForce =
                            curLineContacts[segI][nbI].normalContactForce();
                            
			// if (neiSegI != -1)
                        if (mag(curContactForce) > SMALL)
			{
                            const scalar neiSegParam =
                                curLineContacts[segI][nbI].secondBeamZeta();
                            const label neiLowerSegI =
                                curLineContacts[segI][nbI].secondBeamLowerSegment();
                            const label neiUpperSegI =
                                curLineContacts[segI][nbI].secondBeamUpperSegment();
                        
		            const scalar& curContactDistance =
			        curLineContacts[segI][nbI].delta();
			
			    // const vector curContactForce =
			    //     curLineContacts[segI][nbI].normalContactForce();

			    tensor curContactForceDerivative =
			        curLineContacts[segI][nbI]
			       .normalContactForceDerivative();

			    // const vector curAxialContactForce =
			    //     curLineContacts[segI][nbI]
			    //    .axialContactForce();

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
			    
			    // const label nSegments =
			    //     beam.contact().splines()[nbI].nSegments();

			    const vector dRdZeta =
                                beam.contact().splines()[nbI]
			       .paramFirstDerivative(neiSegI, neiSegParam);
			
			    tensor ImNN = tensor::zero;
			    vector N = curContactForce;
			    if (mag(N) > SMALL)
			    {
                                N /= mag(N);
				ImNN = I - (N*N);
			    }

			    const tensor TT = //tensor::zero;
			        (dRdZeta*dRdZeta)/((dRdZeta & dRdZeta));

			    // Diagonal
			    label globalSegI = start + segI;
			    label ID = 6*globalSegI;

			    tensor FcDn =
			    (
                                mag(curContactForce)*(ImNN - TT)
                               *L[globalSegI]/curContactDistance
			    );

			    // Add frictional component
			    if (mag(curCircumContactForce) > SMALL)
			    {
				curContactForceDerivative +=
			            curCircumContactForceDerivative;

				curContactForceDerivative +=
			            curAxialContactForceDerivative;
			    }

			    //- Contact force derivative (W)
			    A.coeffRef(ID,ID) +=
                                curContactForceDerivative.xx()*L[globalSegI]
			      + FcDn.xx();
			    A.coeffRef(ID,ID+1) +=
			        curContactForceDerivative.xy()*L[globalSegI]
			      + FcDn.xy();
			    A.coeffRef(ID,ID+2) +=
			        curContactForceDerivative.xz()*L[globalSegI]
			      + FcDn.xz();

			    A.coeffRef(ID+1,ID) +=
			        curContactForceDerivative.yx()*L[globalSegI]
			      + FcDn.yx();
			    A.coeffRef(ID+1,ID+1) +=
			        curContactForceDerivative.yy()*L[globalSegI]
			      + FcDn.yy();
			    A.coeffRef(ID+1,ID+2) +=
			        curContactForceDerivative.yz()*L[globalSegI]
			      + FcDn.yz();
			    
			    A.coeffRef(ID+2,ID) +=
			        curContactForceDerivative.zx()*L[globalSegI]
			      + FcDn.zx();
			    A.coeffRef(ID+2,ID+1) +=
			        curContactForceDerivative.zy()*L[globalSegI]
			      + FcDn.zy();
			    A.coeffRef(ID+2,ID+2) +=
			        curContactForceDerivative.zz()*L[globalSegI]
			      + FcDn.zz();

			    //- Contact force derivative (Theta)
			    A.coeffRef(ID,ID+3) +=
                                curCircumContactForceThetaDerivative.xx()*L[globalSegI];
			    A.coeffRef(ID,ID+4) +=
			        curCircumContactForceThetaDerivative.xy()*L[globalSegI];
			    A.coeffRef(ID,ID+5) +=
			        curCircumContactForceThetaDerivative.xz()*L[globalSegI];

			    A.coeffRef(ID+1,ID+3) +=
			        curCircumContactForceThetaDerivative.yx()*L[globalSegI];
			    A.coeffRef(ID+1,ID+4) +=
			        curCircumContactForceThetaDerivative.yy()*L[globalSegI];
			    A.coeffRef(ID+1,ID+5) +=
			        curCircumContactForceThetaDerivative.yz()*L[globalSegI];
			    
			    A.coeffRef(ID+2,ID+3) +=
			        curCircumContactForceThetaDerivative.zx()*L[globalSegI];
			    A.coeffRef(ID+2,ID+4) +=
			        curCircumContactForceThetaDerivative.zy()*L[globalSegI];
			    A.coeffRef(ID+2,ID+5) +=
			        curCircumContactForceThetaDerivative.zz()*L[globalSegI];

			    //- Contact moment derivative
			    A.coeffRef(ID+3,ID+3) +=
                                curContactMomentDerivative.xx()*L[globalSegI];
			    A.coeffRef(ID+3,ID+4) +=
			        curContactMomentDerivative.xy()*L[globalSegI];
			    A.coeffRef(ID+3,ID+5) +=
			        curContactMomentDerivative.xz()*L[globalSegI];

			    A.coeffRef(ID+4,ID+3) +=
			        curContactMomentDerivative.yx()*L[globalSegI];
			    A.coeffRef(ID+4,ID+4) +=
			        curContactMomentDerivative.yy()*L[globalSegI];
			    A.coeffRef(ID+4,ID+5) +=
			        curContactMomentDerivative.yz()*L[globalSegI];
			    
			    A.coeffRef(ID+5,ID+3) +=
			        curContactMomentDerivative.zx()*L[globalSegI];
			    A.coeffRef(ID+5,ID+4) +=
			        curContactMomentDerivative.zy()*L[globalSegI];
			    A.coeffRef(ID+5,ID+5) +=
			        curContactMomentDerivative.zz()*L[globalSegI];

			    
			    // Off-diagonal
			    
			    label neiStart = 0;
			    label i = 0;
			    while(i < nbI)
			    {
			        neiStart +=
				    beam.contact().splines()[i].nSegments();
				i++;
			    }

			    // label globalNeiSegI = neiStart + neiSegI;
                            label globalNeiSeg0 = neiLowerSegI + neiStart;
                            label globalNeiSeg1 = neiUpperSegI + neiStart;
			    scalar w0 = curLineContacts[segI][nbI].secondBeamWeight();
			    scalar w1 = 1.0 - w0;

                            // if (globalNeiSeg0 != globalNeiSeg1)
                            // {
                            //     scalar l =
                            //         0.5*beam.contact().splines()[nbI]
                            //        .segLength(neiLowerSegI)
                            //       + 0.5*beam.contact().splines()[nbI]
                            //        .segLength(neiUpperSegI);

                            //     scalar lx =
                            //         mag(neiSegParam)
			    //        *beam.contact().splines()[nbI]
                            //        .segLength(neiSegI)/2;
                                    
                            //     if (neiSegI == neiLowerSegI)
                            //     {
                            //         w0 = (l-lx)/l;
                            //         w1 = 1.0 - w0;
                            //     }
                            //     else
                            //     {
                            //         w1 = (l-lx)/l;
                            //         w0 = 1.0 - w1;
                            //     }
                            // }
                            
			    // label globalNeiSeg0 = globalNeiSegI;
			    // label globalNeiSeg1 = globalNeiSegI;
			    // scalar w0 = 1.0;
			    // scalar w1 = 0.0;

			    // if
			    // (
                            //     (neiSegI > 0)
			    //  && (neiSegI < (nSegments-1))
			    // )
			    // {
                            //     if (neiSegParam >= 0)
			    //     {
                            //         globalNeiSeg1 = globalNeiSegI + 1;
			    //         scalar l =
			    //             0.5*beam.contact().splines()[nbI]
			    //            .segLength(neiSegI)
			    //           + 0.5*beam.contact().splines()[nbI]
			    //            .segLength(neiSegI+1);
			    //         scalar l0 =
			    //             0.5*neiSegParam
			    //            *beam.contact().splines()[nbI]
			    //            .segLength(neiSegI);
			    //         w0 = (l-l0)/l;
			    //         w1 = 1-w0;
			    //     }
			    //     else
			    //     {
                            //         globalNeiSeg0 = globalNeiSegI - 1;
			    //         scalar l =
			    //             0.5*beam.contact().splines()[nbI]
			    //                .segLength(neiSegI-1)
			    //           + 0.5*beam.contact().splines()[nbI]
			    //                .segLength(neiSegI);
			    //         scalar l1 =
			    //            -0.5*neiSegParam
			    //            *beam.contact().splines()[nbI]
			    //            .segLength(neiSegI);
			    //         w1 = (l-l1)/l;
			    //         w0 = 1-w1;
			    //     }
			    // }

			    label IN0 = 6*globalNeiSeg0;
			    label IN1 = 6*globalNeiSeg1;

			    // IN0
			    
			    //- Contact force derivative (W)
			    A.coeffRef(ID,IN0) -=
			        curContactForceDerivative.xx()
			       *L[globalSegI]*w0 + FcDn.xx()*w0;
			    A.coeffRef(ID,IN0+1) -=
			        curContactForceDerivative.xy()
			       *L[globalSegI]*w0 + FcDn.xy()*w0;
			    A.coeffRef(ID,IN0+2) -=
			        curContactForceDerivative.xz()
			       *L[globalSegI]*w0 + FcDn.xz()*w0;

			    A.coeffRef(ID+1,IN0) -=
			        curContactForceDerivative.yx()
			       *L[globalSegI]*w0 + FcDn.yx()*w0;
			    A.coeffRef(ID+1,IN0+1) -=
			        curContactForceDerivative.yy()
			       *L[globalSegI]*w0 + FcDn.yy()*w0;
			    A.coeffRef(ID+1,IN0+2) -=
			        curContactForceDerivative.yz()
			       *L[globalSegI]*w0 + FcDn.yz()*w0;

			    A.coeffRef(ID+2,IN0) -=
			        curContactForceDerivative.zx()
			       *L[globalSegI]*w0 + FcDn.zx()*w0;
			    A.coeffRef(ID+2,IN0+1) -=
			        curContactForceDerivative.zy()
			       *L[globalSegI]*w0 + FcDn.zy()*w0;
			    A.coeffRef(ID+2,IN0+2) -=
			        curContactForceDerivative.zz()
			       *L[globalSegI]*w0 + FcDn.zz()*w0;

			    //- Contact force derivative (Theta)
			    A.coeffRef(ID,IN0+3) +=
			        curCircumContactForceThetaDerivative.xx()
			       *L[globalSegI]*w0;
			    A.coeffRef(ID,IN0+4) +=
			        curCircumContactForceThetaDerivative.xy()
			       *L[globalSegI]*w0;
			    A.coeffRef(ID,IN0+5) +=
			        curCircumContactForceThetaDerivative.xz()
			       *L[globalSegI]*w0;

			    A.coeffRef(ID+1,IN0+3) +=
			        curCircumContactForceThetaDerivative.yx()
			       *L[globalSegI]*w0;
			    A.coeffRef(ID+1,IN0+4) +=
			        curCircumContactForceThetaDerivative.yy()
			       *L[globalSegI]*w0;
			    A.coeffRef(ID+1,IN0+5) +=
			        curCircumContactForceThetaDerivative.yz()
			       *L[globalSegI]*w0;

			    A.coeffRef(ID+2,IN0+3) +=
			        curCircumContactForceThetaDerivative.zx()
			       *L[globalSegI]*w0;
			    A.coeffRef(ID+2,IN0+4) +=
			        curCircumContactForceThetaDerivative.zy()
			       *L[globalSegI]*w0;
			    A.coeffRef(ID+2,IN0+5) +=
			        curCircumContactForceThetaDerivative.zz()
			       *L[globalSegI]*w0;

			    //- Contact moment derivative
			    A.coeffRef(ID+3,IN0+3) +=
			        curContactMomentDerivative.xx()*L[globalSegI]*w0;
			    A.coeffRef(ID+3,IN0+4) +=
			        curContactMomentDerivative.xy()*L[globalSegI]*w0;
			    A.coeffRef(ID+3,IN0+5) +=
			        curContactMomentDerivative.xz()*L[globalSegI]*w0;

			    A.coeffRef(ID+4,IN0+3) +=
			        curContactMomentDerivative.yx()*L[globalSegI]*w0;
			    A.coeffRef(ID+4,IN0+4) +=
			        curContactMomentDerivative.yy()*L[globalSegI]*w0;
			    A.coeffRef(ID+4,IN0+5) +=
			        curContactMomentDerivative.yz()*L[globalSegI]*w0;

			    A.coeffRef(ID+5,IN0+3) +=
			        curContactMomentDerivative.zx()*L[globalSegI]*w0;
			    A.coeffRef(ID+5,IN0+4) +=
			        curContactMomentDerivative.zy()*L[globalSegI]*w0;
			    A.coeffRef(ID+5,IN0+5) +=
			        curContactMomentDerivative.zz()*L[globalSegI]*w0;
			    
			    // IN1
			    //- Contact force derivative (W)
			    A.coeffRef(ID,IN1) -=
			        curContactForceDerivative.xx()
			       *L[globalSegI]*w1 + FcDn.xx()*w1;
			    A.coeffRef(ID,IN1+1) -=
			        curContactForceDerivative.xy()
			       *L[globalSegI]*w1 + FcDn.xy()*w1;
			    A.coeffRef(ID,IN1+2) -=
			        curContactForceDerivative.xz()
			       *L[globalSegI]*w1 + FcDn.xz()*w1;
			    
			    A.coeffRef(ID+1,IN1) -=
			        curContactForceDerivative.yx()
			       *L[globalSegI]*w1 + FcDn.yx()*w1;
			    A.coeffRef(ID+1,IN1+1) -=
			        curContactForceDerivative.yy()
			       *L[globalSegI]*w1 + FcDn.yy()*w1;
			    A.coeffRef(ID+1,IN1+2) -=
			        curContactForceDerivative.yz()
			       *L[globalSegI]*w1 + FcDn.yz()*w1;
			    
			    A.coeffRef(ID+2,IN1) -=
			        curContactForceDerivative.zx()
			       *L[globalSegI]*w1 + FcDn.zx()*w1;
			    A.coeffRef(ID+2,IN1+1) -=
			        curContactForceDerivative.zy()
			       *L[globalSegI]*w1 + FcDn.zy()*w1;
			    A.coeffRef(ID+2,IN1+2) -=
			        curContactForceDerivative.zz()
			       *L[globalSegI]*w1 + FcDn.zz()*w1;
			    
			    //- Contact force derivative (Theta)
			    A.coeffRef(ID,IN1+3) +=
			        curCircumContactForceThetaDerivative.xx()
			       *L[globalSegI]*w1;
			    A.coeffRef(ID,IN1+4) +=
			        curCircumContactForceThetaDerivative.xy()
			       *L[globalSegI]*w1;
			    A.coeffRef(ID,IN1+5) +=
			        curCircumContactForceThetaDerivative.xz()
			       *L[globalSegI]*w1;

			    A.coeffRef(ID+1,IN1+3) +=
			        curCircumContactForceThetaDerivative.yx()
			       *L[globalSegI]*w1;
			    A.coeffRef(ID+1,IN1+4) +=
			        curCircumContactForceThetaDerivative.yy()
			       *L[globalSegI]*w1;
			    A.coeffRef(ID+1,IN1+5) +=
			        curCircumContactForceThetaDerivative.yz()
			       *L[globalSegI]*w1;

			    A.coeffRef(ID+2,IN1+3) +=
			        curCircumContactForceThetaDerivative.zx()
			       *L[globalSegI]*w1;
			    A.coeffRef(ID+2,IN1+4) +=
			        curCircumContactForceThetaDerivative.zy()
			       *L[globalSegI]*w1;
			    A.coeffRef(ID+2,IN1+5) +=
			        curCircumContactForceThetaDerivative.zz()
			       *L[globalSegI]*w1;
			    
			    //- Contact moment derivative
			    A.coeffRef(ID+3,IN1+3) +=
			        curContactMomentDerivative.xx()*L[globalSegI]*w1;
			    A.coeffRef(ID+3,IN1+4) +=
			        curContactMomentDerivative.xy()*L[globalSegI]*w1;
			    A.coeffRef(ID+3,IN1+5) +=
			        curContactMomentDerivative.xz()*L[globalSegI]*w1;

			    A.coeffRef(ID+4,IN1+3) +=
			        curContactMomentDerivative.yx()*L[globalSegI]*w1;
			    A.coeffRef(ID+4,IN1+4) +=
			        curContactMomentDerivative.yy()*L[globalSegI]*w1;
			    A.coeffRef(ID+4,IN1+5) +=
			        curContactMomentDerivative.yz()*L[globalSegI]*w1;

			    A.coeffRef(ID+5,IN1+3) +=
			        curContactMomentDerivative.zx()*L[globalSegI]*w1;
			    A.coeffRef(ID+5,IN1+4) +=
			        curContactMomentDerivative.zy()*L[globalSegI]*w1;
			    A.coeffRef(ID+5,IN1+5) +=
			        curContactMomentDerivative.zz()*L[globalSegI]*w1;
			}
		    }
		}
	    }

	    start += nSeg;
	}
	}

        // Point contact
        {
            forAll(beam.contact().pointContacts(), pcI)
            {
                labelList bI(2, -1);
                labelList nbI(2, -1);
                labelList segI(2, -1);
                labelList neiSegI(2, -1);
                labelList globalSegI(2, -1);
                labelList globalNeiSegI(2, -1);
                scalarField zeta(2, -2);
                scalarField neiZeta(2, -2);
                scalarField contactDistance(2, 0);
                vectorList contactForce(2, vector::zero);
                tensorList contactForceDerivative(2, tensor::zero);

                vectorList tangentialContactForce(2, vector::zero);
                tensorList tangContactForceDerivative(2, tensor::zero);
                    
                //Owner
                {
                    bI[0] = beam.contact().pointContacts()[pcI].firstBeam();
                    segI[0] =
                        beam.contact().pointContacts()[pcI].firstBeamSegment();
                    zeta[0] =
                        beam.contact().pointContacts()[pcI].firstBeamZeta();
                    label start = 0;
                    label i = 0;
                    while(i < bI[0])
                    {
                        start += beam.contact().splines()[i].nSegments();
                        i++;
                    }
                    globalSegI[0] = start + segI[0];
                    contactDistance[0] =
                        beam.contact().pointContacts()[pcI]
                       .delta();
                    contactForce[0] =
                        beam.contact().pointContacts()[pcI]
                       .normalContactForce();
                    contactForceDerivative[0] =
                        beam.contact().pointContacts()[pcI]
                       .normalContactForceDerivative();

                    tangentialContactForce[0] = 
                        beam.contact().pointContacts()[pcI]
                       .firstBeamTangContactForce();

                    tangContactForceDerivative[0] = 
                        beam.contact().pointContacts()[pcI]
                       .firstBeamTangContactForceDerivative()
                      - beam.contact().pointContacts()[pcI]
                       .secondBeamTangContactForceDerivative();
                }

                // Neighbour
                {
                    bI[1] = beam.contact().pointContacts()[pcI].secondBeam();
                    segI[1] =
                        beam.contact().pointContacts()[pcI].secondBeamSegment();
                    zeta[1] =
                        beam.contact().pointContacts()[pcI].secondBeamZeta();
                    label neiStart = 0;
                    label i = 0;
                    while(i < bI[1])
                    {
                        neiStart += beam.contact().splines()[i].nSegments();
                        i++;
                    }
                    globalSegI[1] = neiStart + segI[1];
                    contactDistance[1] = contactDistance[0];
                    contactForce[1] = -contactForce[0];
                    contactForceDerivative[1] = contactForceDerivative[0];

                    tangentialContactForce[1] = 
                        beam.contact().pointContacts()[pcI]
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
                                mesh.boundary()[beam.endPatchIndex(bI[sI])]
                               .faceCells();

                            label bFaceID =
                                findIndex(faceCells, globalSegI[sI]);

                            DR = zeta[sI]
                               *dRdS.boundaryField()
                                [beam.endPatchIndex(bI[sI])][bFaceID]
                               /dc.boundaryField()
                                [beam.endPatchIndex(bI[sI])][bFaceID];
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
                                mesh.boundary()[beam.startPatchIndex(bI[sI])]
                               .faceCells();

                            label bFaceID =
                                findIndex(faceCells, globalSegI[sI]);

                            DR = zeta[sI]
                               *dRdS.boundaryField()
                                [beam.startPatchIndex(bI[sI])][bFaceID]
                               /dc.boundaryField()
                                [beam.startPatchIndex(bI[sI])][bFaceID];
                        }
                        else
                        {
                            DR = 0.5*zeta[sI]*dRdS.internalField()[faceID]
                               /dc.internalField()[faceID];
                        }
                    }

                    // Diagonal
                    label ID = 6*globalSegI[sI];

                    // label globalSeg0 = globalSegI;
                    label globalSeg1 = globalSegI[sI];
                    scalar w0 = 1;
                    scalar w1 = 0;

                    label nSegments =
                        beam.contact().splines()[bI[sI]].nSegments();
                    scalarField L =
                        beam.contact().splines()[bI[sI]].segLengths();

                    if // this should be extended
                    (
                        (segI[sI] > 0)
                     && (segI[sI] < (nSegments-1))
                    )
                    {
                        if (zeta[sI] >= 0)
                        {
                            globalSeg1 = globalSegI[sI] + 1;
                            scalar l = 0.5*L[segI[sI]] + 0.5*L[segI[sI]+1];
                            scalar l0 = zeta[sI]*L[segI[sI]]/2;
                            w0 = (l-l0)/l;
                            w1 = 1-w0;
                        }
                        else
                        {
                            globalSeg1 = globalSegI[sI] - 1;
                            scalar l = 0.5*L[segI[sI]-1] + 0.5*L[segI[sI]];
                            scalar l1 = -zeta[sI]*L[segI[sI]]/2;
                            w0 = (l-l1)/l;
                            w1 = 1-w0;
                        }
                    }

                    label ID1 = 6*globalSeg1;

                    vector N = contactForce[sI];
                    tensor ImNN = tensor::zero;
                    if (mag(N) > SMALL)
                    {
                        N /= mag(N);
                        ImNN = I - (N*N);
                    }

                    //
                    vector dRdZeta =
                        beam.contact().splines()[bI[sI]]
                       .paramFirstDerivative(segI[sI], zeta[sI]);
                    vector neiDRdZeta =
                        beam.contact().splines()[nbI[sI]]
                       .paramFirstDerivative(neiSegI[sI], neiZeta[sI]);
                    vector dR2dZeta2 =
                        beam.contact().splines()[bI[sI]]
                       .paramSecondDerivative(segI[sI], zeta[sI]);
                    vector neiDR2dZeta2 =
                        beam.contact().splines()[nbI[sI]]
                       .paramSecondDerivative(neiSegI[sI], neiZeta[sI]);

                    vector delta = contactDistance[sI]*N;

                    scalarSquareMatrix M(2, 0.0);
                    M[0][0] = (dRdZeta & neiDRdZeta);
                    M[0][1] = ((delta & neiDR2dZeta2) - (neiDRdZeta & neiDRdZeta));
                    M[1][0] = ((delta & dR2dZeta2) + (dRdZeta & dRdZeta));
                    M[1][1] = -(neiDRdZeta & dRdZeta);

                    scalarSquareMatrix invM = M.LUinvert();

                    tensor TT =
                      - (dRdZeta*(invM[0][0]*neiDRdZeta + invM[0][1]*dRdZeta))
                      + (neiDRdZeta*(invM[1][0]*neiDRdZeta + invM[1][1]*dRdZeta));

                    // TT = tensor::zero;

                    tensor FcDn =
                    (
                        mag(contactForce[sI])
                       *(ImNN + TT)/contactDistance[sI]
                    );
                    // FcDn = tensor::zero;

                    tensor DMCoeff =
                    (
                        spinTensor(DR)
                      & (contactForceDerivative[sI] + FcDn)
                    );
                    // DMCoeff = tensor::zero;

                    //-- w0
                    //- Force
                    A.coeffRef(ID,ID) +=
                        w0*contactForceDerivative[sI].xx()
                      + w0*tangContactForceDerivative[sI].xx()
                      + w0*FcDn.xx();
                    A.coeffRef(ID,ID+1) +=
                        w0*contactForceDerivative[sI].xy()
                      + w0*tangContactForceDerivative[sI].xy()
                      + w0*FcDn.xy();
                    A.coeffRef(ID,ID+2) +=
                        w0*contactForceDerivative[sI].xz()
                      + w0*tangContactForceDerivative[sI].xz()
                      + w0*FcDn.xz();

                    A.coeffRef(ID+1,ID) +=
                        w0*contactForceDerivative[sI].yx()
                      + w0*tangContactForceDerivative[sI].yx()
                      + w0*FcDn.yx();
                    A.coeffRef(ID+1,ID+1) +=
                        w0*contactForceDerivative[sI].yy()
                      + w0*tangContactForceDerivative[sI].yy()
                      + w0*FcDn.yy();
                    A.coeffRef(ID+1,ID+2) +=
                        w0*contactForceDerivative[sI].yz()
                      + w0*tangContactForceDerivative[sI].yz()
                      + w0*FcDn.yz();

                    A.coeffRef(ID+2,ID) +=
                        w0*contactForceDerivative[sI].zx()
                      + w0*tangContactForceDerivative[sI].zx()
                      + w0*FcDn.zx();
                    A.coeffRef(ID+2,ID+1) +=
                        w0*contactForceDerivative[sI].zy()
                      + w0*tangContactForceDerivative[sI].zy()
                      + w0*FcDn.zy();
                    A.coeffRef(ID+2,ID+2) +=
                        w0*contactForceDerivative[sI].zz()
                      + w0*tangContactForceDerivative[sI].zz()
                      + w0*FcDn.zz();

                    //- Moment
                    A.coeffRef(ID+3,ID) +=
                        w0*DMCoeff.xx();
                    A.coeffRef(ID+3,ID+1) +=
                        w0*DMCoeff.xy();
                    A.coeffRef(ID+3,ID+2) +=
                        w0*DMCoeff.xz();

                    A.coeffRef(ID+4,ID) +=
                        w0*DMCoeff.yx();
                    A.coeffRef(ID+4,ID+1) +=
                        w0*DMCoeff.yy();
                    A.coeffRef(ID+4,ID+2) +=
                        w0*DMCoeff.yz();
                    
                    A.coeffRef(ID+5,ID) +=
                        w0*DMCoeff.zx();
                    A.coeffRef(ID+5,ID+1) +=
                        w0*DMCoeff.zy();
                    A.coeffRef(ID+5,ID+2) +=
                        w0*DMCoeff.zz();

                    //-- w1
                    //- Force
                    A.coeffRef(ID,ID1) +=
                        w1*contactForceDerivative[sI].xx()
                      + w1*tangContactForceDerivative[sI].xx()
                      + w1*FcDn.xx();
                    A.coeffRef(ID,ID1+1) +=
                        w1*contactForceDerivative[sI].xy()
                      + w1*tangContactForceDerivative[sI].xy()
                      + w1*FcDn.xy();
                    A.coeffRef(ID,ID1+2) +=
                        w1*contactForceDerivative[sI].xz()
                      + w1*tangContactForceDerivative[sI].xz()
                      + w1*FcDn.xz();

                    A.coeffRef(ID+1,ID1) +=
                        w1*contactForceDerivative[sI].yx()
                      + w1*tangContactForceDerivative[sI].yx()
                      + w1*FcDn.yx();
                    A.coeffRef(ID+1,ID1+1) +=
                        w1*contactForceDerivative[sI].yy()
                      + w1*tangContactForceDerivative[sI].yy()
                      + w1*FcDn.yy();
                    A.coeffRef(ID+1,ID1+2) +=
                        w1*contactForceDerivative[sI].yz()
                      + w1*tangContactForceDerivative[sI].yz()
                      + w1*FcDn.yz();
                    
                    A.coeffRef(ID+2,ID1) +=
                        w1*contactForceDerivative[sI].zx()
                      + w1*tangContactForceDerivative[sI].zx()
                      + w1*FcDn.zx();
                    A.coeffRef(ID+2,ID1+1) +=
                        w1*contactForceDerivative[sI].zy()
                      + w1*tangContactForceDerivative[sI].zy()
                      + w1*FcDn.zy();
                    A.coeffRef(ID+2,ID1+2) +=
                        w1*contactForceDerivative[sI].zz()
                      + w1*tangContactForceDerivative[sI].zz()
                      + w1*FcDn.zz();

                    //- Moment
                    A.coeffRef(ID+3,ID1) +=
                        w1*DMCoeff.xx();
                    A.coeffRef(ID+3,ID1+1) +=
                        w1*DMCoeff.xy();
                    A.coeffRef(ID+3,ID1+2) +=
                        w1*DMCoeff.xz();

                    A.coeffRef(ID+4,ID1) +=
                        w1*DMCoeff.yx();
                    A.coeffRef(ID+4,ID1+1) +=
                        w1*DMCoeff.yy();
                    A.coeffRef(ID+4,ID1+2) +=
                        w1*DMCoeff.yz();

                    A.coeffRef(ID+5,ID1) +=
                        w1*DMCoeff.zx();
                    A.coeffRef(ID+5,ID1+1) +=
                        w1*DMCoeff.zy();
                    A.coeffRef(ID+5,ID1+2) +=
                        w1*DMCoeff.zz();

                    // Off-diagonal
                    label globalNeiSeg0 = globalNeiSegI[sI];
                    label globalNeiSeg1 = globalNeiSegI[sI];
                    w0 = 0.5;
                    w1 = 0.5;

                    label nNeiSegments =
                        beam.contact().splines()[nbI[sI]].nSegments();
                    scalarField neiL =
                        beam.contact().splines()[nbI[sI]].segLengths();

                    if
                    (
                        (neiSegI[sI] > 0)
                     && (neiSegI[sI] < (nNeiSegments-1))
                    )
                    {
                        if (neiZeta[sI] >= 0)
                        {
                            globalNeiSeg1 = globalNeiSegI[sI] + 1;
                            scalar l = 0.5*neiL[neiSegI[sI]]
                              + 0.5*neiL[neiSegI[sI]+1];
                            scalar l0 = neiZeta[sI]*neiL[neiSegI[sI]]/2;
                            w0 = (l-l0)/l;
                            w1 = 1-w0;
                        }
                        else
                        {
                            globalNeiSeg0 = globalNeiSegI[sI] - 1;
                            scalar l = 0.5*neiL[neiSegI[sI]-1]
                              + 0.5*neiL[neiSegI[sI]];
                            scalar l1 = -neiZeta[sI]*neiL[neiSegI[sI]]/2;
                            w1 = (l-l1)/l;
                            w0 = 1-w1;
                        }
                    }

                    label IN0 = 6*globalNeiSeg0;
                    label IN1 = 6*globalNeiSeg1;

                    //-- IN0
                    //- Force
                    A.coeffRef(ID,IN0) -=
                        w0*contactForceDerivative[sI].xx()
                      + w0*tangContactForceDerivative[sI].xx()
                      + w0*FcDn.xx();
                    A.coeffRef(ID,IN0+1) -=
                        w0*contactForceDerivative[sI].xy()
                      + w0*tangContactForceDerivative[sI].xy()
                      + w0*FcDn.xy();
                    A.coeffRef(ID,IN0+2) -=
                        w0*contactForceDerivative[sI].xz()
                      + w0*tangContactForceDerivative[sI].xz()
                      + w0*FcDn.xz();

                    A.coeffRef(ID+1,IN0) -=
                        w0*contactForceDerivative[sI].yx()
                      + w0*tangContactForceDerivative[sI].yx()
                      + w0*FcDn.yx();
                    A.coeffRef(ID+1,IN0+1) -=
                        w0*contactForceDerivative[sI].yy()
                      + w0*tangContactForceDerivative[sI].yy()
                      + w0*FcDn.yy();
                    A.coeffRef(ID+1,IN0+2) -=
                        w0*contactForceDerivative[sI].yz()
                      + w0*tangContactForceDerivative[sI].yz()
                      + w0*FcDn.yz();

                    A.coeffRef(ID+2,IN0) -=
                        w0*contactForceDerivative[sI].zx()
                      + w0*tangContactForceDerivative[sI].zx()
                      + w0*FcDn.zx();
                    A.coeffRef(ID+2,IN0+1) -=
                        w0*contactForceDerivative[sI].zy()
                      + w0*tangContactForceDerivative[sI].zy()
                      + w0*FcDn.zy();
                    A.coeffRef(ID+2,IN0+2) -=
                        w0*contactForceDerivative[sI].zz()
                      + w0*tangContactForceDerivative[sI].zz()
                      + w0*FcDn.zz();

                    //- Moment
                    A.coeffRef(ID+3,IN0) -=
                        w0*DMCoeff.xx();
                    A.coeffRef(ID+3,IN0+1) -=
                        w0*DMCoeff.xy();
                    A.coeffRef(ID+3,IN0+2) -=
                        w0*DMCoeff.xz();

                    A.coeffRef(ID+4,IN0) -=
                        w0*DMCoeff.yx();
                    A.coeffRef(ID+4,IN0+1) -=
                        w0*DMCoeff.yy();
                    A.coeffRef(ID+4,IN0+2) -=
                        w0*DMCoeff.yz();

                    A.coeffRef(ID+5,IN0) -=
                        w0*DMCoeff.zx();
                    A.coeffRef(ID+5,IN0+1) -=
                        w0*DMCoeff.zy();
                    A.coeffRef(ID+5,IN0+2) -=
                        w0*DMCoeff.zz();
                        
                    //-- IN1
                    //- Force
                    A.coeffRef(ID,IN1) -=
                        w1*contactForceDerivative[sI].xx()
                      + w1*tangContactForceDerivative[sI].xx()
                      + w1*FcDn.xx();
                    A.coeffRef(ID,IN1+1) -=
                        w1*contactForceDerivative[sI].xy()
                      + w1*tangContactForceDerivative[sI].xy()
                      + w1*FcDn.xy();
                    A.coeffRef(ID,IN1+2) -=
                        w1*contactForceDerivative[sI].xz()
                      + w1*tangContactForceDerivative[sI].xz()
                      + w1*FcDn.xz();

                    A.coeffRef(ID+1,IN1) -=
                        w1*contactForceDerivative[sI].yx()
                      + w1*tangContactForceDerivative[sI].yx()
                      + w1*FcDn.yx();
                    A.coeffRef(ID+1,IN1+1) -=
                        w1*contactForceDerivative[sI].yy()
                      + w1*tangContactForceDerivative[sI].yy()
                      + w1*FcDn.yy();
                    A.coeffRef(ID+1,IN1+2) -=
                        w1*contactForceDerivative[sI].yz()
                      + w1*tangContactForceDerivative[sI].yz()
                      + w1*FcDn.yz();

                    A.coeffRef(ID+2,IN1) -=
                        w1*contactForceDerivative[sI].zx()
                      + w1*tangContactForceDerivative[sI].zx()
                      + w1*FcDn.zx();
                    A.coeffRef(ID+2,IN1+1) -=
                        w1*contactForceDerivative[sI].zy()
                      + w1*tangContactForceDerivative[sI].zy()
                      + w1*FcDn.zy();
                    A.coeffRef(ID+2,IN1+2) -=
                        w1*contactForceDerivative[sI].zz()
                      + w1*tangContactForceDerivative[sI].zz()
                      + FcDn.zz();
                        
                    //- Moment
                    A.coeffRef(ID+3,IN1) -=
                        w1*DMCoeff.xx();
                    A.coeffRef(ID+3,IN1+1) -=
                        w1*DMCoeff.xy();
                    A.coeffRef(ID+3,IN1+2) -=
                        w1*DMCoeff.xz();

                    A.coeffRef(ID+4,IN1) -=
                        w1*DMCoeff.yx();
                    A.coeffRef(ID+4,IN1+1) -=
                        w1*DMCoeff.yy();
                    A.coeffRef(ID+4,IN1+2) -=
                        w1*DMCoeff.yz();

                    A.coeffRef(ID+5,IN1) -=
                        w1*DMCoeff.zx();
                    A.coeffRef(ID+5,IN1+1) -=
                        w1*DMCoeff.zy();
                    A.coeffRef(ID+5,IN1+2) -=
                        w1*DMCoeff.zz();

                    ////// Tangential force direction corrector
                    if (false)
                    {
                        scalar delta;
                        label ID = 6*globalSegI[sI];
                        label ID1 = 6*globalSeg1;
                        label IND = globalNeiSegI[sI];

                        if (globalSeg1 > globalSegI[sI])
                        {
                            delta = (0.5*L[segI[sI]] + 0.5*L[segI[sI]+1]);                        
                        }
                        else
                        {
                            delta = -(0.5*L[segI[sI]-1] + 0.5*L[segI[sI]]);
                            // Info << "false" << endl;
                        }

                        //
                        A.coeffRef(ID,ID) -= mag(tangentialContactForce[sI])/delta;
                        A.coeffRef(ID+1,ID+1) -= mag(tangentialContactForce[sI])/delta; 
                        A.coeffRef(ID+2,ID+2) -= mag(tangentialContactForce[sI])/delta;
                        
                        A.coeffRef(ID,ID1) += mag(tangentialContactForce[sI])/delta;
                        A.coeffRef(ID+1,ID1+1) += mag(tangentialContactForce[sI])/delta; 
                        A.coeffRef(ID+2,ID1+2) += mag(tangentialContactForce[sI])/delta;

                        //
                        A.coeffRef(IND,ID) += mag(tangentialContactForce[sI])/delta;
                        A.coeffRef(IND+1,ID+1) += mag(tangentialContactForce[sI])/delta; 
                        A.coeffRef(IND+2,ID+2) += mag(tangentialContactForce[sI])/delta;

                        A.coeffRef(IND,ID1) -= mag(tangentialContactForce[sI])/delta;
                        A.coeffRef(IND+1,ID1+1) -= mag(tangentialContactForce[sI])/delta; 
                        A.coeffRef(IND+2,ID1+2) -= mag(tangentialContactForce[sI])/delta;
                    }
                }
            }
        }
    }

    // // Apply conical pulley contact
    // #include "applyConicalPulleyContact.H"
}

// ************************************************************************* //
