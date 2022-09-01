/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
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

Application
    moveBeam

Description
    Transforms the mesh points in the polyMesh directory according to the
    translate, rotate and scale options.

Usage
    Options are:

    -cellZone name

    -translate vector
        Translates the points by the given vector,

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IStringStream.H"
#include "RodriguesRotation.H"
#include "cylindricalCS.H"

using namespace Foam;
using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::validOptions.insert("cellZone", "name");
    argList::validOptions.insert("translate", "vector");
    argList::validOptions.insert("rotateAlongVector", "(vector angleInDegree)");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    if (args.options().empty())
    {
        FatalErrorIn(args.executable())
            << "No options supplied, please use "
            << "-cellZone and -translate options"
            << exit(FatalError);
    }

    fileName meshDir = polyMesh::meshSubDir;

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
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
    
    
    // Translation options

    if (args.optionFound("translate"))
    {
        vector transVector(args.optionLookup("translate")());

        Info<< "Translating beam points by " << transVector << endl;

        const label zoneID = mesh.cellZones().findZoneID(cellZoneName);

        const cellZone& cz = mesh.cellZones()[zoneID];

        const labelListList& pc = mesh.pointCells();

        forAll(pc, pI)
        {
            if (cz.whichCell(pc[pI][0]) > -1)
            {
                points[pI] += transVector;
            }
        }
        
        // points += transVector;
    }
    // else
    // {
    //     FatalErrorIn(args.executable())
    //         << "Option translation is not specified"
    //         << exit(FatalError);
    // }

    // Rotation options

    if (args.optionFound("rotateAlongVector"))
    {
        vector rotationAxis;
        scalar rotationAngle;

        args.optionLookup("rotateAlongVector")()
            >> rotationAxis
            >> rotationAngle;

        tensor T = RodriguesRotation(rotationAxis, rotationAngle);

        Info << "Rotating points by " << T << endl;

        const label zoneID = mesh.cellZones().findZoneID(cellZoneName);

        const cellZone& cz = mesh.cellZones()[zoneID];

        const labelListList& pc = mesh.pointCells();

        forAll(pc, pI)
        {
            if (cz.whichCell(pc[pI][0]) > -1)
            {
                points[pI] = transform(T, points[pI]);
            }
        }
        
        // points = transform(T, points);
    }


    
    Info << "Writing points into directory " << points.path() << nl << endl;
    points.write();

    return 0;
}


// ************************************************************************* //
