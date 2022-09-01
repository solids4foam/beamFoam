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

Description
    Path calculation is performed using procedure described in:
    http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "faceList.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "emptyPolyPatch.H"
#include "Xfer.H"
#include "IOdictionary.H"
#include "cellZone.H"

#include "conicalPulley.H"
#include "OFstream.H"
#include "DynamicField.H"
#include "cylindricalCS.H"
#include "HermiteSpline.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("totalPathLength");

#   include "setRootCase.H"
#   include "createTime.H"

    // Read total path length
    scalar totalPathLength
    (
        readScalar(IStringStream(args.additionalArgs()[0])())
    );

    
    IOdictionary beamProperties
    (
        IOobject
        (
            "beamProperties",
            fileName(runTime.caseConstant()),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE // must be AUTO_WRITE
        )
    );

    PtrList<conicalPulley> conicalPulleys;
    
    // Read conical pulleys
    if (beamProperties.found("conicalPulleys"))
    {
        Info << "Found conical pulleys" << endl;
      
        const PtrList<entry> entries
        (
            beamProperties.lookup("conicalPulleys")
        );

        label nPulleys = entries.size();

        Info << "nPulleys: " << nPulleys << endl;

        if (nPulleys < 1)
        {
            FatalErrorIn(args.executable())
                << "No specified conical pulleys."
                << exit(FatalError);
        }
        
        conicalPulleys.setSize(nPulleys);

        forAll(entries, pulleyI)
        {
            conicalPulleys.set
            (
                pulleyI,
                new conicalPulley(entries[pulleyI].dict())
            );
        }
    }
    else
    {
        FatalErrorIn(args.executable())
            << "No specified conical pulleys."
            << exit(FatalError);
    }

    // Create path
    DynamicField<vector> points;
    DynamicField<vector> tangents;

    //- Add starting point
    {
        vector axis2 =
        (
            conicalPulleys[0].axis()
           ^conicalPulleys[0].polarAxis()
        );

        vector T0 =
            conicalPulleys[0].origin()
          + conicalPulleys[0].maxRadius()*axis2;

        if
        (
            (
                (axis2^conicalPulleys[0].polarAxis())
              & conicalPulleys[0].axis()
            )*conicalPulleys[0].rotDir()
          < 0
        )
        {
            T0 = conicalPulleys[0].origin()
                - conicalPulleys[0].maxRadius()*axis2;
        }

        vector tangent0 = conicalPulleys[0].polarAxis();

        points.append(T0);
        tangents.append(tangent0);
    }

    //- Add remining points
    for (label i=1; i < conicalPulleys.size(); i++)
    {
        if
        (
            mag
            (
                conicalPulleys[i-1].axis()
              & conicalPulleys[i].axis()
            )
          > (1-SMALL)
        )
        {
            coordinateSystem cs
            (
                "cs",
                conicalPulleys[i-1].origin(),
                conicalPulleys[i-1].axis(),
                conicalPulleys[i-1].polarAxis()
            );

            vector C0 =
                cs.localPosition(conicalPulleys[i-1].origin());
            vector C1 =
                cs.localPosition(conicalPulleys[i].origin());

            scalar a = C0.x();
            scalar b = C0.y();

            scalar c = C1.x();
            scalar d = C1.y();

            scalar r0 = conicalPulleys[i-1].maxRadius();
            scalar r1 = conicalPulleys[i].maxRadius();

            scalar xP = (c*r0 + a*r1)/(r0+r1);
            scalar yP = (d*r0 + b*r1)/(r0+r1);

            scalar xt1 =
            (
                sqr(r0)*(xP-a)
              + r0*(yP-b)*::sqrt(sqr(xP-a)+sqr(yP-b)-sqr(r0))
            )
           /(sqr(xP-a)+sqr(yP-b)) + a;

            scalar yt1 =
            (
                sqr(r0)*(yP-b)
              - r0*(xP-a)*::sqrt(sqr(xP-a)+sqr(yP-b)-sqr(r0))
            )
           /(sqr(xP-a)+sqr(yP-b)) + b;

            scalar xt2 =
            (
                sqr(r0)*(xP-a)
              - r0*(yP-b)*::sqrt(sqr(xP-a)+sqr(yP-b)-sqr(r0))
            )
           /(sqr(xP-a)+sqr(yP-b)) + a;

            scalar yt2 =
            (
                sqr(r0)*(yP-b)
              + r0*(xP-a)*::sqrt(sqr(xP-a)+sqr(yP-b)-sqr(r0))
            )
           /(sqr(xP-a)+sqr(yP-b)) + b;

            // scalar s0 = (b-yt1)*(yP-yt1)/(xt1-a)/(xty-xP);

            scalar xt3 =
            (
                sqr(r1)*(xP-c)
              + r1*(yP-d)*::sqrt(sqr(xP-c)+sqr(yP-d)-sqr(r1))
            )
           /(sqr(xP-c)+sqr(yP-d)) + c;

            scalar yt3 =
            (
                sqr(r1)*(yP-d)
              - r1*(xP-c)*::sqrt(sqr(xP-c)+sqr(yP-d)-sqr(r1))
            )
           /(sqr(xP-c)+sqr(yP-d)) + d;

            scalar xt4 =
            (
                sqr(r1)*(xP-c)
              - r1*(yP-d)*::sqrt(sqr(xP-c)+sqr(yP-d)-sqr(r1))
            )
           /(sqr(xP-c)+sqr(yP-d)) + c;

            scalar yt4 =
            (
                sqr(r1)*(yP-d)
              + r1*(xP-c)*::sqrt(sqr(xP-c)+sqr(yP-d)-sqr(r1))
            )
           /(sqr(xP-c)+sqr(yP-d)) + d;

            vector T1(xt1, yt1, 0);
            vector T2(xt2, yt2, 0);

            vector T3(xt3, yt3, 0);
            vector T4(xt4, yt4, 0);
        
            vector d1 = T3 - C1;
            d1 /= mag(d1) + SMALL;

            vector tan1 = T3 - T1;
            tan1 /= mag(tan1) + SMALL;

            vector d2 = T4 - C1;
            d2 /= mag(d2) + SMALL;

            vector tan2 = T4 - T2;
            tan2 /= mag(tan2) + SMALL;

            if
            (
                (
                    (d1^tan1)
                  & cs.localVector
                    (
                        conicalPulleys[i].axis()
                    )
                )*conicalPulleys[i].rotDir()
              > 0
            )
            {
                vector gT1 = cs.globalPosition(T1);
                vector gT3 = cs.globalPosition(T3);
                vector gTan1 = cs.globalVector(tan1);
            
                // First add three points along the arc
                {
                    vector gC0 = conicalPulleys[i-1].origin();

                    label li = points.size()-1;
                
                    vector Tm = (points[li]-gC0) + (gT1-gC0);
                    Tm /= mag(Tm) + SMALL;
                    Tm *= r0;
                
                    // Check rot direction
                    if
                    (
                        (
                            ((points[li] - gC0)^Tm)
                          & conicalPulleys[i-1].axis()
                        )*conicalPulleys[i-1].rotDir()
                      < 0
                    )
                    {
                        Tm *= -1;
                    }
                
                    Tm += gC0;

                    vector tangentm = gT1 - points[li];
                    tangentm /= mag(tangentm) + SMALL;

                    vector Tm0 = (points[li]-gC0) + (Tm-gC0);
                    Tm0 /= mag(Tm0) + SMALL;
                    Tm0 *= r0;
                    Tm0 += gC0;
                    vector tangentm0 = Tm - points[li];
                    tangentm0 /= mag(tangentm0) + SMALL;
                
                    vector Tm1 = (gT1-gC0) + (Tm-gC0);
                    Tm1 /= mag(Tm1) + SMALL;
                    Tm1 *= r0;
                    Tm1 += gC0;
                    vector tangentm1 = gT1 - Tm;
                    tangentm1 /= mag(tangentm1) + SMALL;
                
                    points.append(Tm0);
                    tangents.append(tangentm0);
                
                    points.append(Tm);
                    tangents.append(tangentm);
                
                    points.append(Tm1);
                    tangents.append(tangentm1);
                }

                // Add tangent line points
                points.append(gT1);
                tangents.append(gTan1);
                
                points.append(gT3);
                tangents.append(gTan1);
            }
            else
            {
                vector gT2 = cs.globalPosition(T2);
                vector gT4 = cs.globalPosition(T4);
                vector gTan2 = cs.globalVector(tan2);
            
                // First add three points along the previous arc
                {
                    vector gC0 = conicalPulleys[i-1].origin();

                    label li = points.size()-1;
                
                    vector Tm = (points[li]-gC0) + (gT2-gC0);
                    Tm /= mag(Tm) + SMALL;
                    Tm *= r0;

                    // Check rot direction
                    if
                    (
                        (
                            ((points[li] - gC0)^Tm)
                          & conicalPulleys[i-1].axis()
                        )*conicalPulleys[i-1].rotDir()
                      < 0
                    )
                    {
                        Tm *= -1;
                    }
                
                    Tm += gC0;

                    vector tangentm = gT2 - points[li];
                    tangentm /= mag(tangentm) + SMALL;

                    vector Tm0 = (points[li]-gC0) + (Tm-gC0);
                    Tm0 /= mag(Tm0) + SMALL;
                    Tm0 *= r0;
                    Tm0 += gC0;
                    vector tangentm0 = Tm - points[li];
                    tangentm0 /= mag(tangentm0) + SMALL;
                
                    vector Tm1 = (gT2-gC0) + (Tm-gC0);
                    Tm1 /= mag(Tm1) + SMALL;
                    Tm1 *= r0;
                    Tm1 += gC0;
                    vector tangentm1 = gT2 - Tm;
                    tangentm1 /= mag(tangentm1) + SMALL;

                
                    points.append(Tm0);
                    tangents.append(tangentm0);

                    points.append(Tm);
                    tangents.append(tangentm);
                
                    points.append(Tm1);
                    tangents.append(tangentm1);
                }

                points.append(gT2);
                tangents.append(gTan2);

                points.append(gT4);
                tangents.append(gTan2);
            }
        }
        else // Add transition path segment
        {
            // Add last point of the previous pulley
            {
                label lastPulley = i-1;
        
                vector axis2 =
                (
                    conicalPulleys[lastPulley].axis()
                   ^conicalPulleys[lastPulley].polarAxis()
                );

                vector Te = conicalPulleys[lastPulley].maxRadius()*axis2;

                // Check rot direction
                label lastPoint = points.size()-1;
                if
                (
                    (
                        (
                            (
                                points[lastPoint]
                              - conicalPulleys[lastPulley].origin()
                            ) ^ Te
                        )
                      & conicalPulleys[lastPulley].axis()
                    )*conicalPulleys[lastPulley].rotDir()
                  < 0
                )
                {
                    Te *= -1;
                }

                Te += conicalPulleys[lastPulley].origin();
                vector tangente = conicalPulleys[lastPulley].polarAxis();

                // Add arc central point and tangent
                {
                    vector gC0 = conicalPulleys[lastPulley].origin();
                    scalar r0 = conicalPulleys[lastPulley].maxRadius();
                    label li = points.size()-1;

                    vector Tm = (points[li]-gC0) + (Te-gC0);
                    Tm /= mag(Tm) + SMALL;
                    Tm *= r0;
                    Tm += gC0;
                    
                    vector tangentm = Te - points[li];
                    tangentm /= mag(tangentm);
                    
                    points.append(Tm);
                    tangents.append(tangentm);
                }
                
                points.append(Te);
                tangents.append(tangente);
            }

            // Add straight transition path segment
            {
                vector axis2 =
                (
                    conicalPulleys[i].axis()
                   ^conicalPulleys[i].polarAxis()
                );

                vector T0 =
                    conicalPulleys[i].origin()
                  + conicalPulleys[i].maxRadius()*axis2;

                if
                (
                    (
                        (axis2^conicalPulleys[i].polarAxis())
                      & conicalPulleys[i].axis()
                    )*conicalPulleys[i].rotDir()
                  < 0
                )
                {
                    T0 = conicalPulleys[i].origin()
                      - conicalPulleys[i].maxRadius()*axis2;
                }

                vector tangent0 = conicalPulleys[i].polarAxis();
                
                points.append(T0);
                tangents.append(tangent0);
            }
        }
    }

    //- Add last point
    {
        label lastPulley = conicalPulleys.size()-1;
        
        vector axis2 =
        (
            conicalPulleys[lastPulley].axis()
           ^conicalPulleys[lastPulley].polarAxis()
        );

        vector Te = conicalPulleys[lastPulley].maxRadius()*axis2;

        // Check rot direction
        label lastPoint = points.size()-1;
        if
        (
            (
                ((points[lastPoint] - conicalPulleys[lastPulley].origin())^Te)
              & conicalPulleys[lastPulley].axis()
            )*conicalPulleys[lastPulley].rotDir()
          < 0
        )
        {
            Te *= -1;
        }

        Te += conicalPulleys[lastPulley].origin();
        vector tangente = conicalPulleys[lastPulley].polarAxis();

        // Add arc central point and tangent
        {
            vector gC0 = conicalPulleys[lastPulley].origin();
            scalar r0 = conicalPulleys[lastPulley].maxRadius();
            label li = points.size()-1;

            vector Tm = (points[li]-gC0) + (Te-gC0);
            Tm /= mag(Tm) + SMALL;
            Tm *= r0;
            Tm += gC0;

            vector tangentm = Te - points[li];
            tangentm /= mag(tangentm);

            points.append(Tm);
            tangents.append(tangentm);
        }

        points.append(Te);
        tangents.append(tangente);
    }

    //Add straight part of the path at the outlet
    {
        HermiteSpline currentPath(points, tangents);
        scalar curLength = currentPath.length();
        scalar dl = totalPathLength - curLength;
        if (dl<0)
        {
            FatalErrorIn(args.executable())
                << "Total path length shuld be larger then "
                << curLength << exit(FatalError);            
        }
        
        label lastPulley = conicalPulleys.size()-1;

        vector DL =  dl*conicalPulleys[lastPulley].polarAxis();

        points.append(points[points.size()-1]+DL);
        tangents.append(tangents[tangents.size()-1]);
    }

    // Write path data to file
    {
        fileName conicalPulleysDir(runTime.caseConstant()/"conicalPulleys");
        mkDir(conicalPulleysDir);

        fileName pathFileName("conicalPulleysPath");
        OFstream pathFile(conicalPulleysDir/pathFileName);

        pathFile << points << endl << tangents << endl;
    }

    // Write path data to vtk file
    {
        HermiteSpline spline(points, tangents);

        label n = 1000;
        vectorField pts(n+1);
        scalar L = spline.length();
        scalar DL = L/n;

        Info << "L: " << L << endl;

        for (label i=0; i<=n; i++)
        {
            pts[i] = spline.position(i*DL);
        }

        fileName conicalPulleysDir(runTime.caseConstant()/"conicalPulleys");
        mkDir(conicalPulleysDir);

        fileName vtkFileName("conicalPulleysPath.vtk");
        
        OFstream vtkFile(conicalPulleysDir/vtkFileName);
        
        vtkFile << "# vtk DataFile Version 2.0" << nl
            << vtkFileName << nl
            << "ASCII" << nl
            << "DATASET UNSTRUCTURED_GRID" << nl
            // << "DATASET POLYDATA" << nl
            << "POINTS " << pts.size() << " float" << nl;
        
        forAll(pts, pointI)
        {
            vtkFile << pts[pointI].x() << ' '
                    << pts[pointI].y() << ' '
                    << pts[pointI].z() << nl;
        }

        // Write cells
        label nCells = pts.size()-1;
        vtkFile << "\nCELLS " << nCells << ' ' << nCells*3 << nl;

        for (label i=1; i<pts.size(); i++)
        {
            vtkFile << 2 << ' ' << i-1 << ' ' << i << nl;
        }

        // Write cell types
        vtkFile << "\nCELL_TYPES " << nCells << nl;

        for (label i=1; i<=nCells; i++)
        {
            vtkFile << 3 << nl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
