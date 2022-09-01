/*---------------------------------------------------------------------------* \
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

\*---------------------------------------------------------------------------*/

#include "beamContactModel.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::beamContactModel::updateToroidalPulleyContacts()
{
    label nBeams = beam_.mesh().cellZones().size();
    label nPulleys = beam_.toroidalPulleys().size();

    if (nPulleys)
    {
        // Move pulleys first
        forAll(beam_.toroidalPulleys(), pulleyI)
        {
            toroidalPulley& curPulley =
                const_cast<toroidalPulley&>
                (
                    beam_.toroidalPulleys()[pulleyI]
                );

            curPulley.move(beam_.runTime().timeOutputValue());
        }
        
        // if (beam_.runTime().value() < (beam_.startToRelaxTime()+SMALL))
        {
            // Update nearest points
            for(label bI=0; bI<nBeams; bI++)
            {
                const vectorField& curCentres = splines_[bI].midPoints();

                for (label segI=0; segI<splines_[bI].nSegments(); segI++)
                {    
                    for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                    {
                        vector nearestPoint(vector::zero);

                        nearestPoint =
                            beam_.toroidalPulleys()[pulleyI]
                           .nearestPoint(curCentres[segI]);

                        toroidalPulleyContacts_[bI][segI][pulleyI].set
                        (
                            bI,
                            segI,
                            0,
                            pulleyI,
                            nearestPoint,
                            beam_.runTime().timeIndex()
                        );
                    }
                }
            }
        }
    }
    else
    {
        if (toroidalPulleyContacts_.size())
        {
            toroidalPulleyContacts_.clear();
        }
    }
}

void Foam::beamContactModel::updateToroidalPulleyContactForces()
{
    label nBeams = beam_.mesh().cellZones().size();
    label nPulleys = beam_.toroidalPulleys().size();

    if (nPulleys)
    {
        label start = 0;
        activeToroidalPulleyContact_ = -1.0;
        toroidalPulleyContactForceRatio_ = vector::zero;
        toroidalPulleyNormalContactForce_ = vector::zero;

        for(label bI=0; bI<nBeams; bI++)
        {
            toroidalPulleyContactGaps_[bI] = scalarPair(GREAT, GREAT);

            scalarField avgNormalGap(nPulleys, 0);
            labelList nContactSegments(nPulleys, 0);

            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                label globalSegI = start + segI;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    toroidalPulleyContacts_[bI][segI][pulleyI].updateForce
                    (
                        beam_,
                        splines_,
                        normalModel_(),
                        frictionModel_(),
                        lowerContactAngleLimit_,
                        augmentedLagrangian_
                    );

                    // Store normal gap for each pulley
                    scalar curMinNormalGap =
                        // min
                        (
                            toroidalPulleyContacts_[bI][segI][pulleyI]
                           .normalGap()
                        );

                    if
                    (
                        // max
                        (
                            mag
                            (
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .normalContactForce()
                            )
                        )
                      > SMALL
                    )
                    {
                        scalar curAvgNormalGap =
                            // average
                            (
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .normalGap()
                            );

                        avgNormalGap[pulleyI] += curAvgNormalGap;
                        nContactSegments[pulleyI] += 1;
                    }

                    if
                    (
                        curMinNormalGap
                      < toroidalPulleyContactGaps_[bI][pulleyI].first()
                    )
                    {
                        toroidalPulleyContactGaps_[bI][pulleyI].first() =
                            curMinNormalGap;
                    }

                    if
                    (
                        // max
                        (
                            mag
                            (
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .normalContactForce()
                            )
                        )
                      > SMALL
                    )
                    // if
                    // (
                    //     toroidalPulleyContacts_[bI][segI][pulleyI]
                    //    .activeFriction()
                    // )
                    {
                        // scalar normalContactForce =
                        //     mag
                        //     (
                        //         toroidalPulleyContacts_[bI][segI][pulleyI]
                        //        .normalContactForce()
                        //     );

                        // Info << bI << ", " << segI << ", " << pulleyI
                        //      << ", " << normalContactForce << endl;
                        
                        activeToroidalPulleyContact_[globalSegI] = pulleyI;
                    }

                    // for (label cpI=0; cpI<2; cpI++)
                    {
                        if
                        (
                            toroidalPulleyContacts_[bI][segI][pulleyI]
                           .activeFriction()
                        )
                        {
                            scalar normalContactForce =
                                mag
                                (
                                    toroidalPulleyContacts_[bI][segI][pulleyI]
                                   .normalContactForce()
                                );

                            normalContactForce +=
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .smallFcn_;

                            vector T = splines_[bI].dRdS()[segI];
                            T /= mag(T) + SMALL;
                  
                            scalar frictionalContactForce =
                                // mag
                                (
                                    toroidalPulleyContacts_[bI][segI][pulleyI]
                                   .frictionalContactMoment() & T
                                )/beam_.R(bI);

                            scalar ratio = //frictionalContactForce;
                                frictionalContactForce
                               /normalContactForce;
                            
                            // if (cpI == 0)
                            {
                                toroidalPulleyContactForceRatio_[globalSegI]
                               .x() = ratio;

                                toroidalPulleyNormalContactForce_[globalSegI]
                               .x() = normalContactForce;
                            }
                        }
                        else
                        {
                            vector normalContactForce =
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .normalContactForce();

                            // normalContactForce +=
                            //     toroidalPulleyContacts_[bI][segI][pulleyI]
                            //    .smallFcn_;

                            if (mag(normalContactForce) > SMALL)
                            {
                                toroidalPulleyNormalContactForce_[globalSegI] =
                                    normalContactForce;
                                // vector(normalContactForce, 0, 0);

                                // Info << globalSegI << ", "
                                //      << toroidalPulleyNormalContactForce_[globalSegI] << endl;
                            }
                        }
                    }
                }
            }

            forAll(avgNormalGap, pulleyI)
            {
                reduce(avgNormalGap[pulleyI], sumOp<scalar>());
                reduce(nContactSegments[pulleyI], sumOp<scalar>());
                
                avgNormalGap[pulleyI] /= nContactSegments[pulleyI] + SMALL;
                
                toroidalPulleyContactGaps_[bI][pulleyI].second() =
                    avgNormalGap[pulleyI];

                reduce
                (
                    toroidalPulleyContactGaps_[bI][pulleyI].first(),
                    minOp<scalar>()
                );
            }

            start++;
        }

        // Calculate current load of the pulley bearing
        start = 0;
        toroidalPulleyLoads_ = vector::zero;
        for(label bI=0; bI<nBeams; bI++)
        {
            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    if
                    (
                        toroidalPulleyContacts_[bI][segI][pulleyI]
                       .activeFriction()
                    )
                    {
                        // Calculate load of the pulley bearing
                        vector curContactForce =
                            toroidalPulleyContacts_[bI][segI][pulleyI]
                           .oldNormalContactForce();
                          // + toroidalPulleyContacts_[bI][segI][pulleyI]
                          //  .frictionalContactForce()[0]
                          // + toroidalPulleyContacts_[bI][segI][pulleyI]
                          //  .frictionalContactForce()[1];

                        const vector& curPulleyAxis =
                            beam_.toroidalPulleys()[pulleyI].axis();

                        curContactForce -=
                            curPulleyAxis*(curPulleyAxis&curContactForce);

                        toroidalPulleyLoads_[pulleyI] += curContactForce;
                    }
                }
            }

            start++;
        }

        // Calculate current bearing friction moment
        scalarField Mb(nPulleys, 0);
        forAll(Mb, pulleyI)
        {
            Mb[pulleyI] =
                beam_.toroidalPulleys()[pulleyI].bearingFrictionCoeff()
               *mag(toroidalPulleyLoads_[pulleyI])
               *beam_.toroidalPulleys()[pulleyI].meanBearingDiameter()/2;
        }

        // Calculate distributed axial loading of the beam
        // due to pulley bearing friction        
        start = 0;
        scalarField sumLR(nPulleys, 0);
        for(label bI=0; bI<nBeams; bI++)
        {
            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                label globalSegI = start + segI;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    // for (label cpI=0; cpI<2; cpI++)
                    {
                        scalar curNormalContactForce =
                            mag
                            (
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .oldNormalContactForce()
                            );

                        const vector& pcc =
                            toroidalPulleyContacts_[bI][segI][pulleyI]
                           .pulleyContactCoordinate();
                        
                        scalar smallFcn =
                            toroidalPulleyContacts_[bI][segI][pulleyI].smallFcn_;

                        if (curNormalContactForce > smallFcn)
                        {
                            sumLR[pulleyI] += pcc.x()*beam_.L()[globalSegI];
                        }
                    }
                }
            }

            start++;
        }
        scalarField qa = Mb/(sumLR + SMALL);

        // Update axial frictional contact force
        start = 0;
        for(label bI=0; bI<nBeams; bI++)
        {
            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    // for (label cpI=0; cpI<2; cpI++)
                    {
                        scalar curNormalContactForce =
                            mag
                            (
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .oldNormalContactForce()
                            );

                        scalar smallFcn =
                            toroidalPulleyContacts_[bI][segI][pulleyI].smallFcn_;

                        if (curNormalContactForce > smallFcn)
                        {
                            vector T = splines_[bI].dRdS()[segI];
                            T /= mag(T) + SMALL;

                            const vector& pcc =
                                toroidalPulleyContacts_[bI][segI][pulleyI]
                               .pulleyContactCoordinate();

                            const vector contactPoint =
                                beam_.toroidalPulleys()[pulleyI].position(pcc);

                            const vector& origin =
                                beam_.toroidalPulleys()[pulleyI].origin();

                            const vector R = contactPoint - origin;

                            const vector& axis =
                                beam_.toroidalPulleys()[pulleyI].axis();

                            const label& rotDir =
                                beam_.toroidalPulleys()[pulleyI].rotDir();

                            T *= -sign(rotDir*(axis & (R ^ T)));

                            toroidalPulleyContacts_[bI][segI][pulleyI]
                           .axialFrictionalContactForce() = qa[pulleyI]*T;
                        }
                        else
                        {
                            toroidalPulleyContacts_[bI][segI][pulleyI]
                           .axialFrictionalContactForce() = vector::zero;
                        }
                    }
                }
            }

            start++;
        }        
    }
}

// ************************************************************************* //
