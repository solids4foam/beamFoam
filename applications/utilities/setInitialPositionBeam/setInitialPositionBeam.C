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
    setInitialPositionBeam

Description
    This utility sets the initial configuration of the beams:

    There are two configurations of beam:
    1. Reference Config:
    - Beam mesh when created always lies such that the length direction of beam is aligned to x-axis.
    - The tangents at the face centres of the beam cells are set to (1 0 0) direction.

    2. Initial Configuration of beam: 
    - To start with any other orientation of the beam, one needs to set the position vector of the beam curve.
    - The tangents at the face centres of the beam need to be set from (1 0 0) to the cross-section orientation
        of the beam spatial curve.
    
    * valid arguments are:
    * -cellZone : this command grabs the specififed subset of mesh cells
    * -translate : translate the beam position by the specified vector (x y z)
    * -rotateAlongVector : rotate the beam about a vector with a specified angle - arg reqd ((x y z) angleInDegrees)
    * -rotate : rotate the beam s

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"
#include "boundBox.H"
#include "pointFields.H"

//   newly added
#include "argList.H"
#include "ReadFields.H"
#include "IStringStream.H"
#include "transform.H"


using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    argList::validOptions.insert("cellZone", "name");
    argList::validOptions.insert("translate", "vector");
    argList::validOptions.insert("rotateAlongVector", "(vector angleInDegree)");
    argList::validOptions.insert("rotate", "(vector vector)"); //AT-Added

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    if (args.options().empty())
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use "
            << "-cellZone and -translate options "
			<< "to set the initial position of the beam\n"
			<< "-rotateAlongVector option is optional to use"
            << exit(FatalError);
    }

    surfaceVectorField  refTangent
    (
        IOobject
        (
            "refTangent",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
         ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );


    // Read cell zone
    word cellZoneName;

    if (args.optionFound("cellZone"))
    {
        cellZoneName = args.optionRead<word>("cellZone");

        Info << "Beam cell zone: " << cellZoneName << endl;
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Option cellZone is not specified"
            << exit(FatalError);
    }

    // Translation and rotation options
    if (args.optionFound("translate"))
    {
        vector transVector(args.optionLookup("translate")());

	    tensor T = tensor::I;

        if (args.optionFound("rotate"))
        {
            Pair<vector> n1n2(args.optionLookup("rotate")());
            n1n2[0] /= mag(n1n2[0]);
            n1n2[1] /= mag(n1n2[1]);
            T = rotationTensor(n1n2[0], n1n2[1]);

            Info<< "Rotating points by " << T << endl;
        }

        const label zoneID = mesh.cellZones().findZoneID(cellZoneName);


	  Info << "zoneID" << zoneID << endl;

	 if (zoneID == -1)
	{
	    FatalErrorIn
	    (
		"setInitialPositionBeam application utility"
	    )
		<< "Problem in beam cellZone"
		<< "\nzoneID of beam: " << cellZoneName << " is " << zoneID
		<< "\nProvide the beam name without the hyphen in front "
		<< "e.g. beam_0"
		<< abort(FatalError);
	}

	vectorField& refWfI = refWf.primitiveFieldRef();
	vectorField& refWI = refW.primitiveFieldRef();
	vectorField& refTangentI = refTangent.primitiveFieldRef();
	tensorField& refLambdaI = refLambda.primitiveFieldRef();
	tensorField& refRMI = refRM.primitiveFieldRef();
	const vectorField& CfI = mesh.Cf().primitiveField();
	const vectorField& CI = mesh.C().primitiveField();

	vectorField newCfI ((T & CfI));
	const labelList& nei = mesh.neighbour();

	//- Set all the internal face fields
	forAll(refWfI, faceI)
	{
	    label I = mesh.cellZones().whichZone(nei[faceI]);

	    if (I == zoneID)
	    {
		refWfI[faceI] =  newCfI[faceI] + transVector - CfI[faceI];

		refLambdaI[faceI] = T;

                // PC: refTang shoudl be wrt to the reference frame, not
                // incremental; otherwise, refTang would be inconsistent with
                // refWf and refW
		//refTangentI[faceI] = (refLambdaI[faceI] & refTangentI[faceI]);
		refTangentI[faceI] = (refLambdaI[faceI] & vector(1, 0, 0));
	    }
	}

	vectorField newCI ((T & CI));

	//- Set all the cell-centre reference fields
	forAll(refWI, cellI)
	{
	    label I  = mesh.cellZones().whichZone(cellI);
	    Info << "loop I=" << I << endl;

	    if (I == zoneID)
	    {
		refWI[cellI]  = newCI[cellI] + transVector - CI[cellI];

		refRMI[cellI] = T;
	    }
	}

	//- Set the reference fields at the boundaries
	forAll(refWf.boundaryField(), patchI)
	{
	    vectorField& pRefWf = refWf.boundaryFieldRef()[patchI];
	    vectorField& pRefW = refW.boundaryFieldRef()[patchI];
	    vectorField& pRefTangent = refTangent.boundaryFieldRef()[patchI];
	    tensorField& pRefLambda = refLambda.boundaryFieldRef()[patchI];
	    tensorField& pRefRM = refRM.boundaryFieldRef()[patchI];

	    const vectorField& pCf = mesh.Cf().boundaryField()[patchI];

	    const vectorField newpCf ((T & pCf));

	    const labelList faceCells =
		mesh.boundary()[patchI].faceCells();

	    forAll(pRefWf, faceI)
	    {
		label I = mesh.cellZones().whichZone(faceCells[faceI]);

		if (I == zoneID)
		{
		    pRefWf[faceI] =  newpCf[faceI] + transVector - pCf[faceI];
		    pRefW[faceI] = pRefWf[faceI];

		    pRefLambda[faceI] = T;
		    pRefRM[faceI] = pRefLambda[faceI];

                    // PC: see comment above
		    // pRefTangent[faceI] = (pRefLambda[faceI] & pRefTangent[faceI]);
		    pRefTangent[faceI] = (pRefLambda[faceI] & vector(1, 0, 0));
		}
	    }
	}
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Option translate is not specified"
            << exit(FatalError);
    }

    refWf.write();
    refW.write();
    refLambda.write();
    refTangent.write();
    refRM.write();

#include "updateMeshPoints.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
