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
    Sets the initial configuration of the beams
    * valid arguments are:
    * -cellZone
    * -translate
    * -rotateAlongVector

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
    argList::validOptions.insert("rotate", "(vector vector)"); //AT-Added- the type of vector needs to be fixed

#   include "setRootCase.H"
#   include "createTime.H"
//#   include "createMesh.H"

	const word regionName = args.optionRead<word>("region");
	Info<<regionName<<endl;
	Foam::fvMesh mesh
(
    Foam::IOobject
    (
        //Foam::fvMesh::defaultRegion,
        regionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);
	//meshPtr=&(time_.lookupObject<fvMesh>(regionName));
	//const fvMesh& mesh= *meshPtr;


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

        // Info<< "Translating beam points by " << transVector << endl;

	if (args.optionFound("rotateAlongVector"))// AT: this needs a minor fix, sth not working here :|||||
	{
	    vector rotationAxis;
	    scalar rotationAngle;

	    args.optionLookup("rotateAlongVector")()
		>> rotationAxis
		>> rotationAngle;
	    scalar rotationAngleRadians = rotationAngle * (3.14159265359/180);
	    //Info << rotationAngleRadians << endl;

	    T = Ra(rotationAxis, rotationAngleRadians);
	    //Info << T << endl;
	    Info << "Rotation tensor\n" << T << endl;
	}
	else
	{
	    T  = tensor::I;
	}
	if (args.optionFound("rotate"))
    {
        Pair<vector> n1n2(args.optionLookup("rotate")());
        n1n2[0] /= mag(n1n2[0]);
        n1n2[1] /= mag(n1n2[1]);
        T = rotationTensor(n1n2[0], n1n2[1]);

        Info<< "Rotating points by " << T << endl;

        //points = transform(T, points);

        //if (args.optionFound("rotateFields"))
        //{
          //  rotateFields(args, runTime, T);
        //}
    }

        const label zoneID = mesh.cellZones().findZoneID(cellZoneName);

		//const label regionName = args.optionRead<word>("region");


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

        //   const cellZone& cz = mesh.cellZones()[zoneID];

        //   const labelListList& pc = mesh.pointCells();

        //   forAll(pc, pI)
        //   {
            //   if (cz.whichCell(pc[pI][0]) > -1)
            //   {
                //   points[pI] += transVector;
            //   }
        //   }

        // points += transVector;

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

		refTangentI[faceI] = (refLambdaI[faceI] & refTangentI[faceI]);
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

		    pRefTangent[faceI] = (pRefLambda[faceI] & pRefTangent[faceI]);
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
    //Info<< refW<<endl;
    refLambda.write();
    //Info<<"some Space"<<endl;
    //Info<< refLambda<<endl;
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
