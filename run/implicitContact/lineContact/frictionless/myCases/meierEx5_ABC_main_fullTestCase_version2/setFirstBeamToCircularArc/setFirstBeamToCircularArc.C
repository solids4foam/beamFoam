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

Application
    setFirstBeamToCircularArc

Description
    Deforming straight beam geometry to the 180deg arc shape

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"
#include "boundBox.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    // Beam cross-section rotation
    /*
    volVectorField  refTheta
    (
        IOobject
        (
            "refTheta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    );
    */

    surfaceVectorField  refTangent
    (
        IOobject
        (
            "refTangent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("x", dimless, vector(1, 0, 0))
    );
    
    // Beam mean line displacement
    surfaceVectorField refWf
    (
        IOobject
        (
            "refWf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );


    // Beam mean line displacement
    volVectorField refW
    (
        IOobject
        (
            "refW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );
    
    // Beam cross-section rotation
    surfaceTensorField refLambda
    (
        IOobject
        (
            "refLambda",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );

    volTensorField refRM
    (
        IOobject
        (
            "refRM",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );
    tensorField& refRMI = refRM.internalField();

    //argList::validOptions.insert("arcAngleInDegrees", "scalar");
    
    //scalar arcAngle = (readScalar(args.optionLookup("arcAngleInDegrees")()))*M_PI/180;
    
    scalar arcAngle = M_PI;
    
    // Calc displacement and rotation
    //boundBox box(mesh.points());
    scalar R0 = 1.0;
    scalar L = arcAngle*R0;
    // Info << "Beam length: " << L << endl;
    
    //Info << "Arc radius: " << R0 << endl;
    
    // vectorField& refThetaI = refTheta.internalField();
    vectorField& refWfI = refWf.internalField();
    vectorField& refWI = refW.internalField();
    tensorField& refLambdaI = refLambda.internalField();
    vectorField& refTangentI = refTangent.internalField();
    const vectorField& CfI = mesh.Cf().internalField();
    const vectorField& CI = mesh.C().internalField();
    
    const labelList& nei = mesh.neighbour();
    
    forAll(refWfI, faceI)
    {
	label I =  mesh.cellZones().whichZone(nei[faceI]);
	
	//Set only first Beam
	if (I == 0)
	{
	    scalar x = CfI[faceI].x();
	    scalar beta = arcAngle*x/L;
	    scalar alpha = 0.5*(2*M_PI - arcAngle);
	    scalar phi = M_PI - (alpha + beta);

	    refWfI[faceI] =
		vector
		(
		    R0*(1 - (::cos(beta) + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::sin(alpha),
		    R0*(1 - (::cos(beta) + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::cos(alpha) + R0*::sin(beta)/::sin(alpha),
		    0
		)
	      - vector(x, 0, 0);

	    refLambdaI[faceI] =
		tensor
		(
		    ::cos(phi), -::sin(phi), 0,
		    ::sin(phi),  ::cos(phi), 0,
		    0,             0,            1
		);

	    refTangentI[faceI] = (refLambdaI[faceI] & refTangentI[faceI]);
	}
    }
    
    
    forAll(CI, cellI)
    {
	label I =  mesh.cellZones().whichZone(nei[cellI]);
	
	//Set only first Beam
	if (I == 0)
	{
	    scalar x = CI[cellI].x();
	    scalar beta = arcAngle*x/L;
	    scalar alpha = 0.5*(2*M_PI - arcAngle);
	    scalar phi = M_PI - (alpha + beta);

	    refRMI[cellI] =
		tensor
		(
		    ::cos(phi), -::sin(phi), 0,
		    ::sin(phi),  ::cos(phi), 0,
		    0,             0,            1
		);
	    
	    refWI[cellI] =
		vector
		(
		    R0*(1 - (::cos(beta) 
		    + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::sin(alpha),
		    R0*(1 - (::cos(beta) 
		    + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::cos(alpha)
		    + R0*::sin(beta)/::sin(alpha),
		    0
		)
	      - vector(x, 0, 0);
       }
    }

    /*
    forAll(refThetaI, cellI)
    {
        scalar x = CI[cellI].x();
        scalar beta = arcAngle*x/L;

        refThetaI[cellI] = vector(0, 0, beta);
    }
    */
    
    forAll(refWf.boundaryField(), patchI)
    {
        // vectorField& pRefTheta = refTheta.boundaryField()[patchI];
        vectorField& pRefWf = refWf.boundaryField()[patchI];
        tensorField& pRefLambda = refLambda.boundaryField()[patchI];
        vectorField& pRefTangent = refTangent.boundaryField()[patchI];
        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];
        
        tensorField& pRefRM = refRM.boundaryField()[patchI];
	
	const labelList faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(pRefWf, faceI)
        {
	    label I =  mesh.cellZones().whichZone(faceCells[faceI]);
	
	    //Set only first Beam
	    if (I == 0)
	    {
		scalar x = pCf[faceI].x();
		scalar beta = arcAngle*x/L;
		scalar alpha = 0.5*(2*M_PI - arcAngle);
		scalar phi = M_PI - (alpha + beta);

		pRefWf[faceI] =
		    vector
		(
		    R0*(1 - (::cos(beta) + 
		    ::sin(beta)*::cos(alpha)/::sin(alpha)))*::sin(alpha),
		    R0*(1 - (::cos(beta) 
		    + ::sin(beta)*::cos(alpha)/::sin(alpha)))*::cos(alpha) 
		    + R0*::sin(beta)/::sin(alpha),
		    0
		)
		  - vector(x, 0, 0);

		pRefLambda[faceI] =
		    tensor
		    (
			::cos(phi), -::sin(phi), 0,
			::sin(phi),  ::cos(phi), 0,
			0,             0,        1
		    );

		// pRefTheta[faceI] = vector(0, 0, beta);

		pRefTangent[faceI] = (pRefLambda[faceI] & pRefTangent[faceI]);
		
		pRefRM[faceI] = pRefLambda[faceI];
	    }
	}
    }
	
    refWf.write();
    refW.write();
    refLambda.write();
    refTangent.write();
    
    refRM.write();

#   include "updateMeshPoints.H"
  
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
