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

void Foam::beamContactModel::updateConicalPulleyContacts()
{
    label nBeams = beam_.mesh().cellZones().size();
    label nPulleys = beam_.conicalPulleys().size();

    if (nPulleys)
    {
        // Move pulleys first
        forAll(beam_.conicalPulleys(), pulleyI)
        {
            conicalPulley& curPulley =
                const_cast<conicalPulley&>
                (
                    beam_.conicalPulleys()[pulleyI]
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
                        vectorField nearestPoints(2, vector::zero);

                        nearestPoints[0] =
                            beam_.conicalPulleys()[pulleyI]
                           .nearestPointNeg(curCentres[segI]);
                        nearestPoints[1] =
                            beam_.conicalPulleys()[pulleyI]
                           .nearestPointPos(curCentres[segI]);

                        conicalPulleyContacts_[bI][segI][pulleyI].set
                        (
                            bI,
                            segI,
                            0,
                            pulleyI,
                            nearestPoints,
                            beam_.runTime().timeIndex()
                        );
                    }
                }
            }
        }
    }
    else
    {
        if (conicalPulleyContacts_.size())
        {
            conicalPulleyContacts_.clear();
        }
    }
}

void Foam::beamContactModel::updateConicalPulleyContactForces()
{
    label nBeams = beam_.mesh().cellZones().size();
    label nPulleys = beam_.conicalPulleys().size();

    if (nPulleys)
    {
        label start = 0;
        activeConicalPulleyContact_ = -1.0;
        conicalPulleyContactForceRatio_ = vector::zero;
        conicalPulleyNormalContactForce_ = vector::zero;
        conicalPulleyAxialContactForce_ = vector::zero;

        for(label bI=0; bI<nBeams; bI++)
        {
            conicalPulleyContactGaps_[bI] = scalarPair(GREAT, GREAT);

            scalarField avgNormalGap(nPulleys, 0);
            labelList nContactSegments(nPulleys, 0);

            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                // label globalSegI = start + segI;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    conicalPulleyContacts_[bI][segI][pulleyI]
                   .updateNormalForce
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
                        min
                        (
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .normalGap()
                        );

                    if
                    (
                        max
                        (
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .normalContactForce()
                            )
                        )
                      > SMALL
                    )
                    {
                        scalar curAvgNormalGap =
                            average
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .normalGap()
                            );

                        avgNormalGap[pulleyI] += curAvgNormalGap;
                        nContactSegments[pulleyI] += 1;
                    }

                    if
                    (
                        curMinNormalGap
                      < conicalPulleyContactGaps_[bI][pulleyI].first()
                    )
                    {
                        conicalPulleyContactGaps_[bI][pulleyI].first() =
                            curMinNormalGap;
                    }
                }
            }

            forAll(avgNormalGap, pulleyI)
            {
                reduce(avgNormalGap[pulleyI], sumOp<scalar>());
                reduce(nContactSegments[pulleyI], sumOp<scalar>());
                
                avgNormalGap[pulleyI] /= nContactSegments[pulleyI] + SMALL;
                
                conicalPulleyContactGaps_[bI][pulleyI].second() =
                    avgNormalGap[pulleyI];

                reduce
                (
                    conicalPulleyContactGaps_[bI][pulleyI].first(),
                    minOp<scalar>()
                );
            }

            start++;
        }


        // Update axial and transversal tangential gap
        start = 0;
        for(label bI=0; bI<nBeams; bI++)
        {
            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                // label globalSegI = start + segI;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    conicalPulleyContacts_[bI][segI][pulleyI]
                   .updateTangentialGap
                    (
                        beam_,
                        splines_,
                        normalModel_(),
                        frictionModel_(),
                        lowerContactAngleLimit_,
                        augmentedLagrangian_
                    );                    
                }
            }

            start++;
        }


        // Update axial and transversal tangential gap
        // for steady state cases (convection)
        if (true)
        {
        start = 0;
        for(label bI=0; bI<nBeams; bI++)
        {
            if (beam_.U(bI) > SMALL)
            {
                // Info << "Forward convection" << endl;
                // Forward convection
                for (label segI=1; segI<splines_[bI].nSegments(); segI++)
                {
                    for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                    {
                        const scalarField& ttgiUpwind =
                            conicalPulleyContacts_[bI][segI-1][pulleyI]
                           .transversalTangentialGapIncrement();

                        scalarField& ttgiCurrent =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .transversalTangentialGapIncrement();

                        const scalarField& atgiUpwind =
                            conicalPulleyContacts_[bI][segI-1][pulleyI]
                           .axialTangentialGapIncrement();

                        scalarField& atgiCurrent =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .axialTangentialGapIncrement();
                        
                        const vectorField& fcnUpwind =
                            conicalPulleyContacts_[bI][segI-1][pulleyI]
                           .oldNormalContactForce();
                        const vectorField& fcnCurrent =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .oldNormalContactForce();

                        scalar smallFcn =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .smallFcn_;

                        for (label cpI=0; cpI<2; cpI++)
                        {
                            label fupwind = 0;
                            label fcurrent = 0;

                            if (mag(fcnUpwind[cpI]) > smallFcn)
                            {
                                fupwind = 1;
                            }

                            if (mag(fcnCurrent[cpI]) > smallFcn)
                            {
                                fcurrent = 1;
                            }

                            ttgiCurrent[cpI] +=
                                fcurrent*(fupwind*ttgiUpwind[cpI]);
                            
                            atgiCurrent[cpI] +=
                                fcurrent*(fupwind*atgiUpwind[cpI]);
                        }
                    }
                }
            }
            else if (beam_.U(bI) < -SMALL)
            {
                // Info << "Backward convection" << endl;
                
                // Backward convection
                for (label segI=(splines_[bI].nSegments()-2); segI>=0; segI--)
                {
                    for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                    {
                        const scalarField& ttgiUpwind =
                            conicalPulleyContacts_[bI][segI+1][pulleyI]
                           .transversalTangentialGapIncrement();
                        
                        scalarField& ttgiCurrent =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .transversalTangentialGapIncrement();
                        
                        const scalarField& atgiUpwind =
                            conicalPulleyContacts_[bI][segI+1][pulleyI]
                           .axialTangentialGapIncrement();
                        
                        scalarField& atgiCurrent =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .axialTangentialGapIncrement();
                        
                        const vectorField& fcnUpwind =
                            conicalPulleyContacts_[bI][segI+1][pulleyI]
                           .oldNormalContactForce();
                        const vectorField& fcnCurrent =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .oldNormalContactForce();
                        
                        scalar smallFcn =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .smallFcn_;

                        for (label cpI=0; cpI<2; cpI++)
                        {
                            label fupwind = 0;
                            label fcurrent = 1;

                            if (mag(fcnUpwind[cpI]) > smallFcn)
                            {
                                fupwind = 1;
                            }

                            if
                            (
                                (mag(fcnCurrent[0]) > smallFcn)
                             && (mag(fcnCurrent[1]) > smallFcn)
                            )
                            {
                                fcurrent = 0;
                            }
                            else if (mag(fcnCurrent[cpI]) < smallFcn)
                            {
                                fupwind = 0;
                            }

                            ttgiCurrent[cpI] +=
                                fcurrent*(fupwind*ttgiUpwind[cpI]);
                            
                            atgiCurrent[cpI] +=
                                fcurrent*(fupwind*atgiUpwind[cpI]);
                        }
                    }
                }
            }

            start++;
        }
        }
        
        // Update axial and transversal contact force and moment
        for(label bI=0; bI<nBeams; bI++)
        {
            conicalPulleyContactGaps_[bI] = scalarPair(GREAT, GREAT);

            scalarField avgNormalGap(nPulleys, 0);
            labelList nContactSegments(nPulleys, 0);

            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                label globalSegI = start + segI;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    conicalPulleyContacts_[bI][segI][pulleyI]
                   .updateFrictionForceAndMoment
                    (
                        beam_,
                        splines_,
                        normalModel_(),
                        frictionModel_(),
                        lowerContactAngleLimit_,
                        augmentedLagrangian_
                    );
                    
                    if
                    (
                        conicalPulleyContacts_[bI][segI][pulleyI]
                       .activeFriction()[0]
                     || conicalPulleyContacts_[bI][segI][pulleyI]
                       .activeFriction()[1]
                    )
                    {
                        activeConicalPulleyContact_[globalSegI] = pulleyI;
                    }

                    for (label cpI=0; cpI<2; cpI++)
                    {
                        scalar normalContactForce =
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .normalContactForce()[cpI]
                            );

                        scalar smallFcn =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .smallFcn_;

                        if (normalContactForce > smallFcn)
                        {
                            if (cpI == 0)
                            {
                                conicalPulleyNormalContactForce_[globalSegI]
                               .x() = normalContactForce;
                            }
                            else
                            {
                                conicalPulleyNormalContactForce_[globalSegI]
                               .y() = normalContactForce;
                            }
                        }

                        if
                        (
                            normalContactForce > smallFcn
                           // true
                           //  conicalPulleyContacts_[bI][segI][pulleyI]
                           // .activeFriction()[cpI]
                        )
                        {                            
                            vector T = splines_[bI].dRdS()[segI];
                            T /= mag(T) + SMALL;

                            scalar frictionalContactForce =
                                // mag
                                (
                                    conicalPulleyContacts_[bI][segI][pulleyI]
                                   .transversalFrictionalContactMoment()[cpI] & T
                                )/beam_.R(bI);

                            scalar axialContactForce =
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .axialFrictionalContactForce()[cpI] & T
                            )/(normalContactForce + smallFcn);

                            scalar ratio = //frictionalContactForce;
                                frictionalContactForce
                               /(normalContactForce + smallFcn);

                            if (cpI == 0)
                            {
                                conicalPulleyContactForceRatio_[globalSegI]
                               .x() = ratio;

                                conicalPulleyAxialContactForce_[globalSegI]
                               .x() = axialContactForce;
                            }
                            else
                            {
                                conicalPulleyContactForceRatio_[globalSegI]
                               .y() = ratio;

                                conicalPulleyAxialContactForce_[globalSegI]
                               .y() = axialContactForce;
                                
                                bool unloading =
                                    conicalPulleyContacts_[bI][segI][pulleyI]
                                   .unloading()[cpI];

                                if (unloading)
                                {
                                    conicalPulleyContactForceRatio_[globalSegI]
                                        .z() = 1.0;
                                }
                                else
                                {
                                    conicalPulleyContactForceRatio_[globalSegI]
                                        .z() = 0.0;
                                }
                            }
                        }
                    }
                }
            }

            start++;
        }

        
        // Calculate current load of the pulley bearing
        start = 0;
        conicalPulleyLoads_ = vector::zero;
        for(label bI=0; bI<nBeams; bI++)
        {
            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                // label globalSegI = start + segI;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    if
                    (
                        conicalPulleyContacts_[bI][segI][pulleyI]
                       .activeFriction()[0]
                     || conicalPulleyContacts_[bI][segI][pulleyI]
                       .activeFriction()[1]
                    )
                    {
                        // Calculate load of the pulley bearing
                        vector curContactForce =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .oldNormalContactForce()[0]
                          + conicalPulleyContacts_[bI][segI][pulleyI]
                           .oldNormalContactForce()[1];
                          // + conicalPulleyContacts_[bI][segI][pulleyI]
                          //  .frictionalContactForce()[0]
                          // + conicalPulleyContacts_[bI][segI][pulleyI]
                          //  .frictionalContactForce()[1];

                        const vector& curPulleyAxis =
                            beam_.conicalPulleys()[pulleyI].axis();

                        curContactForce -=
                            curPulleyAxis*(curPulleyAxis&curContactForce);

                        conicalPulleyLoads_[pulleyI] += curContactForce;
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
                beam_.conicalPulleys()[pulleyI].bearingFrictionCoeff()
               *mag(conicalPulleyLoads_[pulleyI])
               *beam_.conicalPulleys()[pulleyI].meanBearingDiameter()/2;
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
                    for (label cpI=0; cpI<2; cpI++)
                    {
                        scalar curNormalContactForce =
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .oldNormalContactForce()[cpI]
                            );

                        const vector& pcc =
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .pulleyContactCoordinates()[cpI];
                        
                        scalar smallFcn =
                            conicalPulleyContacts_[bI][segI][pulleyI].smallFcn_;

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
                // label globalSegI = start + segI;

                for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                {
                    for (label cpI=0; cpI<2; cpI++)
                    {
                        scalar curNormalContactForce =
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .oldNormalContactForce()[cpI]
                            );

                        scalar smallFcn =
                            conicalPulleyContacts_[bI][segI][pulleyI].smallFcn_;

                        if (curNormalContactForce > smallFcn)
                        {
                            vector T = splines_[bI].dRdS()[segI];
                            T /= mag(T) + SMALL;

                            const vector& pcc =
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .pulleyContactCoordinates()[cpI];

                            const vector contactPoint =
                                beam_.conicalPulleys()[pulleyI].position(pcc);

                            const vector& origin =
                                beam_.conicalPulleys()[pulleyI].origin();

                            const vector R = contactPoint - origin;

                            const vector& axis =
                                beam_.conicalPulleys()[pulleyI].axis();

                            const label& rotDir =
                                beam_.conicalPulleys()[pulleyI].rotDir();

                            T *= -sign(rotDir*(axis & (R ^ T)));

                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .axialBearingFrictionForce()[cpI] = qa[pulleyI]*T;
                        }
                        else
                        {
                            conicalPulleyContacts_[bI][segI][pulleyI]
                           .axialBearingFrictionForce()[cpI] = vector::zero;
                        }
                    }
                }
            }

            start++;
        }
    }
}

// ************************************************************************* //
