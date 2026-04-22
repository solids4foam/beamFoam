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
    makeCircularArc

Description
    Sets the initial configuration of the beam to a circular arc
    of specified radius and arc angle in the global xy-plane only.
    * valid arguments are:
    * -cellZone (word - beam type)
    * -arcRadius (scalar)
    * -arcAngleInDegrees (scalar)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"
#include "boundBox.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("cellZone", "name");
    argList::validOptions.insert("arcRadius", "scalar");
    argList::validOptions.insert("arcAngleInDegrees", "scalar");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (args.options().empty())
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use "
            << "-cellZone, -arcRadius and "
            << "-arcAngleInDegrees options to "
            << "set the circular arc \n"
            << exit(FatalError);
    }

    IOdictionary beamProperties
    (
        IOobject
        (
            "beamProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    surfaceVectorField refTangent
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

    surfaceTensorField refLambdaf
    (
        IOobject
        (
            "refLambdaf",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, tensor::I)
    );

    volTensorField refLambda
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

    word cellZoneName;

    if (args.found("cellZone"))
    {
        cellZoneName = args.get<word>("cellZone");

        Info<< "Beam cell zone: " << cellZoneName << endl;
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Option cellZone is not specified"
            << exit(FatalError);
    }

    const scalar arcAngle =
        readScalar(args.lookup("arcAngleInDegrees")())*M_PI/180;

    const scalar R0 = readScalar(args.lookup("arcRadius")());
    const scalar L = mag(arcAngle)*R0;

    const label zoneID = mesh.cellZones().findZoneID(cellZoneName);

    if (zoneID == -1)
    {
        FatalErrorIn("makeCircularArc application utility")
            << "Problem in beam cellZone"
            << "\nzoneID of beam: " << cellZoneName << " is " << zoneID
            << "\nProvide the beam name without the hyphen in front "
            << "e.g. beam_0"
            << abort(FatalError);
    }

    const PtrList<entry> entries(beamProperties.lookup("beams"));
    DynamicList<label> beamEntryIDs;

    forAll(entries, entryI)
    {
        const string entryName(entries[entryI].keyword());
        label nEntryBeams = 1;

        for (const char c : entryName)
        {
            if (c == '|')
            {
                ++nEntryBeams;
            }
        }

        for (label beamI = 0; beamI < nEntryBeams; ++beamI)
        {
            beamEntryIDs.append(entryI);
        }
    }

    if (zoneID >= beamEntryIDs.size())
    {
        FatalErrorIn("makeCircularArc application utility")
            << "Cannot map cellZone " << cellZoneName
            << " with zoneID " << zoneID
            << " to an entry in constant/beamProperties beams."
            << abort(FatalError);
    }

    const scalar beamLength =
        readScalar(entries[beamEntryIDs[zoneID]].dict().lookup("length"));

    if (mag(beamLength - L) > 0.01)
    {
        FatalErrorIn("makeCircularArc application utility")
            << "Inconsistent length of beam provided with the "
            << "length required to create the arc of the "
            << "specified radius and arcAngle "
            << "\nbeam length provided: " << beamLength
            << "\nbeam length needed: " << L
            << "\nBeam length should be equal to arcRadius*arcAngle"
            << abort(FatalError);
    }

    vectorField& refWfI = refWf.primitiveFieldRef();
    vectorField& refWI = refW.primitiveFieldRef();
    tensorField& refLambdafI = refLambdaf.primitiveFieldRef();
    tensorField& refLambdaI = refLambda.primitiveFieldRef();
    vectorField& refTangentI = refTangent.primitiveFieldRef();

    const vectorField& CfI = mesh.Cf().primitiveField();
    const vectorField& CI = mesh.C().primitiveField();
    const labelList& nei = mesh.neighbour();

    forAll(refWfI, faceI)
    {
        const label I = mesh.cellZones().whichZone(nei[faceI]);

        if (I == zoneID)
        {
            const scalar x = CfI[faceI].x();
            const scalar beta = arcAngle*x/L;
            const scalar alpha = 0.5*(2*M_PI - arcAngle);
            const scalar phi = M_PI - (alpha + beta);

            refWfI[faceI] =
                vector
                (
                    R0
                   *(
                        1
                      - (
                            ::cos(beta)
                          + ::sin(beta)*::cos(alpha)/::sin(alpha)
                        )
                    )
                   *::sin(alpha),
                    R0
                   *(
                        1
                      - (
                            ::cos(beta)
                          + ::sin(beta)*::cos(alpha)/::sin(alpha)
                        )
                    )
                   *::cos(alpha)
                  + R0*::sin(beta)/::sin(alpha),
                    0
                )
              - vector(x, 0, 0);

            refLambdafI[faceI] =
                tensor
                (
                    ::cos(phi), -::sin(phi), 0,
                    ::sin(phi),  ::cos(phi), 0,
                    0,             0,        1
                );

            refTangentI[faceI] = (refLambdafI[faceI] & vector(1, 0, 0));
        }
    }

    forAll(refWI, cellI)
    {
        const label I = mesh.cellZones().whichZone(cellI);

        if (I == zoneID)
        {
            const scalar x = CI[cellI].x();
            const scalar beta = arcAngle*x/L;
            const scalar alpha = 0.5*(2*M_PI - arcAngle);
            const scalar phi = M_PI - (alpha + beta);

            refLambdaI[cellI] =
                tensor
                (
                    ::cos(phi), -::sin(phi), 0,
                    ::sin(phi),  ::cos(phi), 0,
                    0,             0,        1
                );

            refWI[cellI] =
                vector
                (
                    R0
                   *(
                        1
                      - (
                            ::cos(beta)
                          + ::sin(beta)*::cos(alpha)/::sin(alpha)
                        )
                    )
                   *::sin(alpha),
                    R0
                   *(
                        1
                      - (
                            ::cos(beta)
                          + ::sin(beta)*::cos(alpha)/::sin(alpha)
                        )
                    )
                   *::cos(alpha)
                  + R0*::sin(beta)/::sin(alpha),
                    0
                )
              - vector(x, 0, 0);
        }
    }

    forAll(refWf.boundaryField(), patchI)
    {
        vectorField& pRefWf = refWf.boundaryFieldRef()[patchI];
        vectorField& pRefW = refW.boundaryFieldRef()[patchI];
        tensorField& pRefLambdaf = refLambdaf.boundaryFieldRef()[patchI];
        tensorField& pRefLambda = refLambda.boundaryFieldRef()[patchI];
        vectorField& pRefTangent = refTangent.boundaryFieldRef()[patchI];

        const vectorField& pCf = mesh.Cf().boundaryField()[patchI];
        const labelList faceCells = mesh.boundary()[patchI].faceCells();

        forAll(pRefWf, faceI)
        {
            const label I = mesh.cellZones().whichZone(faceCells[faceI]);

            if (I == zoneID)
            {
                const scalar x = pCf[faceI].x();
                const scalar beta = arcAngle*x/L;
                const scalar alpha = 0.5*(2*M_PI - arcAngle);
                const scalar phi = M_PI - (alpha + beta);

                pRefWf[faceI] =
                    vector
                    (
                        R0
                       *(
                            1
                          - (
                                ::cos(beta)
                              + ::sin(beta)*::cos(alpha)/::sin(alpha)
                            )
                        )
                       *::sin(alpha),
                        R0
                       *(
                            1
                          - (
                                ::cos(beta)
                              + ::sin(beta)*::cos(alpha)/::sin(alpha)
                            )
                        )
                       *::cos(alpha)
                      + R0*::sin(beta)/::sin(alpha),
                        0
                    )
                  - vector(x, 0, 0);

                pRefLambdaf[faceI] =
                    tensor
                    (
                        ::cos(phi), -::sin(phi), 0,
                        ::sin(phi),  ::cos(phi), 0,
                        0,             0,        1
                    );

                pRefTangent[faceI] =
                    (pRefLambdaf[faceI] & vector(1, 0, 0));

                pRefLambda[faceI] = pRefLambdaf[faceI];
                pRefW[faceI] = pRefWf[faceI];
            }
        }
    }

    refWf.write();
    refW.write();
    refLambdaf.write();
    refTangent.write();
    refLambda.write();

    #include "updateMeshPoints.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
