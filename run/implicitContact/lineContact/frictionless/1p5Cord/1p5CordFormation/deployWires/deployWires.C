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
    deployWires

Description
    Deploy 5 wires around core wire

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoDPointCorrector.H"
#include "boundBox.H"
#include "pseudoVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    // Calc displacement and rotation
    boundBox box(mesh.points());
    scalar L = box.max().x() - box.min().x();
    Info << "Beam length: " << L << endl;

    // Radius of enveloping cylinder
    scalar R0 = 0.41;
    Info << "Enveloping cylinder radius: " << R0 << endl;

    label nBeams = mesh.pointZones().size();

    pointIOField newPoints
    (
        IOobject
        (
            "points",
            runTime.constant(),
            // runTime.findInstance(polyMesh::meshSubDir, "points"),
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    for (label i=1; i<nBeams; i++)
    {
        const labelList& curBeamPoints = mesh.pointZones()[i];

        forAll(curBeamPoints, pI)
        {
            newPoints[curBeamPoints[pI]] += vector(0, R0, 0);
        }
    }
    
    for (label i=2; i<nBeams; i++)
    {
        const labelList& curBeamPoints = mesh.pointZones()[i];

        scalar DTheta = (360/5)*M_PI/180;
        
        scalar Theta = (i-1)*DTheta;
        
        // Rotation tensor
        tensor A
        (
            1, 0, 0,
            0, ::cos(Theta), -::sin(Theta),
            0, ::sin(Theta),  ::cos(Theta)
        );
        
        forAll(curBeamPoints, pI)
        {
            newPoints[curBeamPoints[pI]] =
                (A & newPoints[curBeamPoints[pI]]);
        }            
    }
        
    Info << "Writing new mesh points ...";
    newPoints.write();
    Info << " done" << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
