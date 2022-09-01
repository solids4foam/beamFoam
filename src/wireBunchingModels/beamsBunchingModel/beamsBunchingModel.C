/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield160
         | foam-extend: Open Source CFD
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

\*---------------------------------------------------------------------------*/

#include "beamsBunchingModel.H"
#include "OStringStream.H"
#include "cubicSpline.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(beamsBunchingModel, 0);
//     defineRunTimeSelectionTable(beamsBunchingModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::beamsBunchingModel::updateContact()
{
    PtrList<cubicSpline> splines(nearestPoints_.size());
    forAll(splines, sI)
    {
        splines.set
        (
            sI,
            new cubicSpline
            (
                beams_[sI].currentBeamPoints(),
                cubicSpline::CLAMPED_CALC,
                vector(0, 0, 0),
                cubicSpline::CLAMPED_CALC,
                vector(0, 0, 0)
            )
        );
    }

    // Update nearest points
    for(label bI=0; bI<nearestPoints_.size(); bI++)
    {
        labelScalarListList& curNearestPoints = nearestPoints_[bI];
        const vectorField& curCentres = splines[bI].midPoints();

        for (label nbI=0; nbI<nearestPoints_.size(); nbI++)
        {
            if (nbI != bI) // No self-contact
            {
                const vectorField& neiPoints = splines[nbI].points();
                const vectorField& neidRdS = splines[nbI].dRdS();

                for (label segI=0; segI<curNearestPoints[nbI].size(); segI++)
                {
                    // Find nearest segment
                    bool found = false;
                    for
                    (
                        label nbSegI=0;
                        nbSegI<beams_[nbI].mesh().nCells();
                        nbSegI++
                    )
                    {
                        scalar testValue =
                        (
                            neidRdS[nbSegI]
                          & (curCentres[segI] - neiPoints[nbSegI])
                        )
                       *(
                            neidRdS[nbSegI+1]
                          & (curCentres[segI] - neiPoints[nbSegI+1])
                        );

                        if (testValue <= 0)
                        {
                            nearestPoints_[bI][nbI][segI] = 
                                splines[nbI].nearestPoint
                                (
                                    nbSegI,
                                    curCentres[segI]
                                );
                            // nearestPoints_[bI][nbI][segI].second() = 0.5;
                            // Info << segI << ", "
                            //      << nearestPoints_[bI][nbI][segI]
                            //      <<  endl;
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        Info << "Can not find nearest segment for"
                             << " seg: " << segI << " of beam: " << bI << endl;
                    }
                }
            }
        }
    }

    // Calculate contact force
    scalar maxResidual = 0;
    for(label bI=0; bI<beams_.size(); bI++)
    {
        const vectorField& C = splines[bI].midPoints();

        const vectorField DW = 
            beams_[bI].currentDisplacementIncrement();
        const tensorField DLambda =
            beams_[bI].currentRotationIncrement(); 

        scalar minGap = GREAT;

        forAll(nearestPoints_[bI], nbI)
        {
            const labelScalarList& curNearestPoints = nearestPoints_[bI][nbI];

            vectorList& curContactForces = contactForces_[bI][nbI];
            vectorList& curFrictionalContactForces =
                frictionalContactForces_[bI][nbI];
            vectorList& curFrictionalContactMoments =
                frictionalContactMoments_[bI][nbI];

            scalarList& curContactGap = contactGap_[bI][nbI];
            // scalarList& curContactOffset = contactOffset_[bI][nbI];

            // if (bI < nbI) // Slave (and no self-contact)
            if (bI != nbI) // No self-contact
            {
                // scalar maxFc = max(mag(curContactForces)) + SMALL;

                const vectorField neiDW = 
                    beams_[nbI].currentDisplacementIncrement();
                const tensorField neiDLambda =
                    beams_[nbI].currentRotationIncrement();

                forAll(curNearestPoints, segI)
                {
                    vector curPoint =
                        splines[nbI].position
                        (
                            curNearestPoints[segI].first(),
                            curNearestPoints[segI].second()
                        );

                    scalar dist = mag(C[segI] - curPoint);

                    vector n = C[segI] - curPoint;
                    n /= mag(n) + SMALL;

                    scalar gap = dist - beams_[bI].R()
                      - beams_[nbI].R();

                    // gap += 0.0010001; // offset

                    if (gap < minGap)
                    {
                        minGap = gap;
                    }
                    
                    scalar curResidual =
                        mag(gap - curContactGap[segI])
                       /beams_[bI].R();
                    
                    curContactGap[segI] = gap;

                    vector fc =
                        contactForce
                        (
                            curContactGap[segI]
                          // - curContactOffset[segI]
                        )*n;

                    // if (curContactGap[segI] > 0)
                    // {
                    //   Info << curContactGap[segI] << ", "
                    //        << fc << endl;
                    // }
                    
                    
                    // scalar curResidual =
                    //     mag(fc-curContactForces[segI])/maxFc;

                    if (curResidual > maxResidual)
                    {
                        maxResidual = curResidual;
                    }

                    curContactForces[segI] =
                        forceRelaxFac_*fc
                      + (1.0 - forceRelaxFac_)*curContactForces[segI];

                    // Caculate friction contact forces
                    if (frictionCoeff_ > SMALL)
                    {
                        if (mag(curContactForces[segI]) > SMALL)
                        {
                            label neiSegI = curNearestPoints[segI].first();
                            scalar neiSegParam =
                                curNearestPoints[segI].second();

                            // Contact point displacement
                            vector DR = -n*beams_[bI].R();
                            tensor curDLambda =
                                0.5*(DLambda[segI]+DLambda[segI+1]);
                            vector curDWt = DW[segI]
                              + ((I-curDLambda.T()) & DR);
                            curDWt -= n*(n & curDWt);

                            // Contact point displacement (neighbour beam)
                            vector neiDR = n*beams_[nbI].R();
                            tensor curNeiDLambda = neiDLambda[neiSegI]
                              + neiSegParam
                               *(neiDLambda[neiSegI+1]-neiDLambda[neiSegI]);
                            vector curNeiDWt = neiDW[neiSegI]
                              + ((I-curNeiDLambda.T()) & neiDR);
                            curNeiDWt -= n*(n & curDWt);

                            vector curSlip = (curDWt - curNeiDWt);
                            curSlip /= mag(curSlip) + SMALL;

                            vector ffc =
                               -frictionCoeff_
                               *mag(curContactForces[segI])*curSlip;
                            
                            curFrictionalContactForces[segI] =
                                forceRelaxFac_*ffc
                              + (1.0 - forceRelaxFac_)
                               *curFrictionalContactForces[segI];


                            curFrictionalContactMoments[segI] =
                                (DR ^ ffc);
                            
                            // curFrictionalContactMoments[segI] =
                            //     forceRelaxFac_
                            //    *(DR ^ curFrictionalContactForces[segI])
                            //   + (1.0 - forceRelaxFac_)
                            //    *curFrictionalContactMoments[segI];
                        }
                    }
                }
            }
            // else if (bI > nbI) // Master (transfering force from slave)
            // {
            //     forAll(curNearestPoints, segI)
            //     {
            //         label neiSegI = curNearestPoints[segI].first();
            //         scalar neiParam = curNearestPoints[segI].second();

            //         const vectorField& curCentres = splines[nbI].midPoints();
            //         vector curPoint = splines[nbI].position(neiSegI, neiParam);

            //         const vectorList& neiContactForces = 
            //             contactForces_[nbI][bI];

            //         if (neiParam <= 0.5)
            //         {
            //             if (neiSegI>0)
            //             {
            //                 vector f0 = neiContactForces[neiSegI-1];
            //                 vector f1 = neiContactForces[neiSegI];

            //                 scalar s01 = mag(curPoint-curCentres[neiSegI-1])
            //                   + mag(curCentres[neiSegI]-curPoint);
            //                 scalar s = mag(curPoint-curCentres[neiSegI-1]);

            //                 curContactForces[segI] = -(f0 + (f1-f0)*s/s01);
            //             }
            //             else
            //             {
            //                 curContactForces[segI] = -neiContactForces[neiSegI];
            //             }
            //         }
            //         else
            //         {
            //             if (neiSegI<(neiContactForces.size()-1))
            //             {
            //                 vector f0 = neiContactForces[neiSegI];
            //                 vector f1 = neiContactForces[neiSegI+1];

            //                 scalar s01 = mag(curPoint-curCentres[neiSegI])
            //                   + mag(curCentres[neiSegI+1]-curPoint);
            //                 scalar s = mag(curPoint-curCentres[neiSegI]);

            //                 curContactForces[segI] = -(f0 + (f1-f0)*s/s01);
            //             }
            //             else
            //             {
            //                 curContactForces[segI] = -neiContactForces[neiSegI];
            //             }
            //         }    
            //     }
            // }
        }

        Info << "Beam: " << bI << ", minGap: " << minGap
             << ", fc: " << contactForce(minGap) << endl;
        // beams_[bI].q().internalField() = q;
    }

    // Info << contactForces_ << endl;
        
    return maxResidual;
}

void Foam::beamsBunchingModel::setContactForces()
{
    for (label bI=0; bI<beams_.size(); bI++)
    {
        const vectorListList& curContactForces = contactForces_[bI];
        const vectorListList& curFrictionalContactForces =
            frictionalContactForces_[bI];
        const vectorListList& curFrictionalContactMoments =
            frictionalContactMoments_[bI];
        
        vectorField& curq = beams_[bI].q().internalField();
        curq = vector::zero;

        vectorField& curm = beams_[bI].m().internalField();
        curm = vector::zero;
        
        forAll(curContactForces, nbI)
        {
            if (nbI!=bI)
            {
                curq += curContactForces[nbI];
                curq += curFrictionalContactForces[nbI];
                curm += curFrictionalContactMoments[nbI];

                Info << "avg Fc: "
                    << average(mag(curContactForces[nbI])) << endl;
            }
        }
    }
}

Foam::scalar Foam::beamsBunchingModel::contactForce(const scalar g) const
{
    scalar fc = 0;

    dimensionedScalar epsilon(lookup("epsilon"));

    dimensionedScalar gBar(lookup("gBar"));
    gBar.value() += SMALL;

    // scalar gap = g;
    scalar gap = g;
    // scalar gap = g + gBar.value();
    
    scalar fBar = epsilon.value()*gBar.value()/2;

    if (gap <= 0)
    {
        fc = fBar - epsilon.value()*gap;
    }
    else if (gap <= gBar.value())
    {
        fc = (epsilon.value()*gBar.value() - fBar)*sqr(gap/gBar.value())
          - epsilon.value()*gap + fBar;
    }
    
    return fc;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::beamsBunchingModel::beamsBunchingModel(Time& runTime)
:
    IOdictionary
    (
        IOobject
        (
            "beamsBunchingProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(runTime),
    beams_(),
    nearestPoints_(),
    contactForces_(),
    frictionalContactForces_(),
    frictionalContactMoments_(),
    contactGap_(),
    contactOffset_(),
    frictionCoeff_(lookupOrDefault<scalar>("frictionCoeff", 0)),
    forceRelaxFac_(readScalar(lookup("forceRelaxationFactor"))),
    gapRelaxFac_(readScalar(lookup("gapRelaxationFactor"))),
    convergenceTol_(readScalar(lookup("convergenceTol"))),
    nCorr_(readInt(lookup("nCorrectors"))),
    residualFilePtr_(NULL),
    active_(lookupOrDefault<bool>("active", true))
{
    label nBeams(readInt(lookup("nBeams")));

    beams_.setSize(nBeams);
    forAll(beams_, bI)
    {
        OStringStream BeamName;
        BeamName() << "beam-" << bI;

        beams_.set
        (
            bI,
            beamModel::New(runTime, word(BeamName.str()))
        );
    }

    // Initialize nearest points
    nearestPoints_.setSize(nBeams);
    forAll(nearestPoints_, bI)
    {
        nearestPoints_.set
        (
            bI,
            new labelScalarListList
            (
                nBeams
            )
        );
        forAll(nearestPoints_[bI], nbI)
        {
            if (nbI != bI)
            {
                nearestPoints_[bI][nbI] =
                    labelScalarList
                    (
                        beams_[bI].mesh().nCells(),
                        labelScalar(-1, 0)
                    );
            }
        }
    }

    // Initialize contact forces
    contactForces_.setSize(nBeams);
    frictionalContactForces_.setSize(nBeams);
    frictionalContactMoments_.setSize(nBeams);
    forAll(contactForces_, bI)
    {
        contactForces_.set
        (
            bI,
            new vectorListList
            (
                nBeams
            )
        );
        frictionalContactForces_.set
        (
            bI,
            new vectorListList
            (
                nBeams
            )
        );
        frictionalContactMoments_.set
        (
            bI,
            new vectorListList
            (
                nBeams
            )
        );
        forAll(contactForces_[bI], nbI)
        {
            if (nbI != bI)
            {
                contactForces_[bI][nbI] =
                    vectorList
                    (
                        beams_[bI].mesh().nCells(),
                        vector::zero
                    );
                frictionalContactForces_[bI][nbI] =
                    vectorList
                    (
                        beams_[bI].mesh().nCells(),
                        vector::zero
                    );
                frictionalContactMoments_[bI][nbI] =
                    vectorList
                    (
                        beams_[bI].mesh().nCells(),
                        vector::zero
                    );
            }
        }
    }

    // Initialize contact gap
    contactGap_.setSize(nBeams);
    forAll(contactGap_, bI)
    {
        contactGap_.set
        (
            bI,
            new scalarListList
            (
                nBeams
            )
        );
        forAll(contactGap_[bI], nbI)
        {
            if (nbI != bI)
            {
                contactGap_[bI][nbI] =
                    scalarList
                    (
                        beams_[bI].mesh().nCells(),
                        GREAT
                    );
            }
        }
    }

    // Initialize contact offset
    contactOffset_.setSize(nBeams);
    forAll(contactOffset_, bI)
    {
        contactOffset_.set
        (
            bI,
            new scalarListList
            (
                nBeams
            )
        );
        forAll(contactOffset_[bI], nbI)
        {
            if (nbI != bI)
            {
                contactOffset_[bI][nbI] =
                    scalarList
                    (
                        beams_[bI].mesh().nCells(),
                        0.0
                    );
            }
        }
    }
    
    // If requested, create the residual file
    if (lookupOrDefault<Switch>("residualFile", false))
    {
        if (Pstream::master())
        {
            Info<< "Creating residual.dat" << endl;
            residualFilePtr_.set
            (
                new OFstream(runTime.path()/"residual.dat")
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::beamsBunchingModel::~beamsBunchingModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::beamsBunchingModel::evolve()
{
    nCorr_ = label(readInt(lookup("nCorrectors")));
    forceRelaxFac_ = scalar(readScalar(lookup("forceRelaxationFactor")));
    convergenceTol_ = scalar(readScalar(lookup("convergenceTol")));

    bool updatePenaltyLawOffset
    (
        lookupOrDefault<bool>
        (
            "updatePenaltyLawOffset",
            false
        )
    );
    
    label iCorr = 0;
    scalar residual = 0;
    do
    {
        if (active_)
        {
            setContactForces();
        }
        
        residual = 0;
        forAll(beams_, bI)
        {
            scalar curRes = beams_[bI].evolve();

            if (curRes > residual)
            {
                residual = curRes;
            }
        }

        // updateContact();
        residual = updateContact();

        Info << "Current contact residual: " << residual << endl;
        Info << "Current force relaxation: " << forceRelaxFac_ << endl;
    }
    while(++iCorr<nCorr_ && residual>convergenceTol_);

    if (residualFilePtr_.valid())
    {
        residualFilePtr_() << runTime().value() << " "
                           << residual << " " << iCorr << endl;
    }

    // Update offset
    if (updatePenaltyLawOffset)
    {
        for (label bI=0; bI<beams_.size(); bI++)
        {
            for (label nbI=0; nbI<beams_.size(); nbI++)
            {
                vectorList& curContactForces = contactForces_[bI][nbI];
                scalarList& curContactGap = contactGap_[bI][nbI];
                scalarList& curContactOffset = contactOffset_[bI][nbI];

                if (bI != nbI)
                {
                    scalarField contact = pos(mag(curContactForces)-SMALL);
                    scalar avgGap =
                       sum(contact*curContactGap)/(sum(contact) + SMALL);

                    Info << average(curContactOffset) << ", " << avgGap << endl;
                
                    if (avgGap < 0)
                    {
                      curContactOffset = curContactOffset - avgGap;
                    }

                }
            }
        }
    }
    
    Info << "iCorr: " << iCorr << ", forceRelaxFactor: "
         << forceRelaxFac_ << endl;
}

void Foam::beamsBunchingModel::updateTotalFields()
{
    forAll(beams_, bI)
    {
        beams_[bI].updateTotalFields();
    }
}

void Foam::beamsBunchingModel::writeFields()
{
    if (runTime().outputTime())
    {
        runTime().write();

        Switch writeVTK = false;
        if (found("writeVTK"))
        {
            writeVTK = Switch(lookup("writeVTK"));
        }

        if (writeVTK)
        {
            forAll(beams_, bI)
            {
                beams_[bI].writeVTK();
            }
        }

        // Write contact force and contact gap for first beam
        for (label bI=0; bI<beams_.size(); bI++)
        {
            // label bI = 0;
            cubicSpline spline
            (
                beams_[bI].points(),
                cubicSpline::CLAMPED_CALC,
                vector(0, 0, 0),
                cubicSpline::CLAMPED_CALC,
                vector(0, 0, 0)
            );
            
            scalarField t = spline.midPointParameters();

            // Info << "First beam length: " << sum(spline.segLengths()) << endl;

            Info << "Total forces for beam " << bI << endl;
            for (label i=0; i<beams_.size(); i++)
            {
                if (i != bI)
                {
                    Info << ' ' << sum(contactForces_[bI][i]);
                }
                else
                {
                    Info << ' ' << 0;
                }
            }
            Info << endl;


            
            
            {
                OStringStream FileName;
                FileName() << "beam-" << bI << "_force.txt";
        
                // OFstream forceFile(runTime().timePath()/"beam-0_force.txt");
                OFstream forceFile
                (
                    runTime().timePath()/word(FileName.str())
                );
                forAll(t, segI)
                {
                    forceFile << t[segI];
                    for (label i=0; i<beams_.size(); i++)
                    {
                        if (i != bI)
                        {
                            forceFile << ' ' << mag(contactForces_[bI][i][segI]);
                        }
                        else
                        {
                            forceFile << ' ' << 0;
                        }
                    }
                    forceFile << endl;
                }
            }

            {
                OStringStream FileName;
                FileName() << "beam-" << bI << "_gap.txt";
        
                OFstream gapFile
                (
                    runTime().timePath()/word(FileName.str())
                );
                
                // OFstream gapFile(runTime().timePath()/"beam-0_gap.txt");
                forAll(t, segI)
                {
                    gapFile << t[segI];
                    for (label i=0; i<beams_.size(); i++)
                    {
                        if (i != bI)
                        {
                            gapFile << ' ' << contactGap_[bI][i][segI];
                        }
                        else
                        {
                            gapFile << ' ' << 0;
                        }
                    }
                    gapFile << endl;
                }
            }
        }
    }
}

// ************************************************************************* //
