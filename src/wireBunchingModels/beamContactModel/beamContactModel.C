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
#include "OStringStream.H"
#include "cubicSpline.H"
#include "OFstream.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(beamContactModel, 0);
//     defineRunTimeSelectionTable(beamContactModel, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::beamContactModel::initializeLineContact
(
    const label bI,
    const label segI,
    const label nbI
)
{
	
    // Info << "initializeLineContact: " <<  bI << ", " << segI << endl;
  
    const vectorField& neiPoints = splines_[nbI].points();
    const vectorField& neidRdS = splines_[nbI].dRdS();
    const vectorField& curCentres = splines_[bI].midPoints();

    labelScalar curNearestPoint(-1, 0);

    // Find nearest segment; Contact Search for the nearest segment
    // in the neighouring beam
    for
    (
        label nbSegI = 0; nbSegI < splines_[nbI].nSegments(); nbSegI++
    )
    {
        scalar testValue =
		(
			neidRdS[nbSegI]
			& (curCentres[segI] - neiPoints[nbSegI])
		)
       *(
			neidRdS[nbSegI + 1]
			& (curCentres[segI] - neiPoints[nbSegI + 1])
		);
		
		// Info << "Contact search testValue " << testValue << endl;

		if (testValue <= 0)
		{
			curNearestPoint = 
				splines_[nbI].nearestPoint
				(
				nbSegI,
				curCentres[segI]
				);
			break;
		}
    }

    label nbSegI = curNearestPoint.first();
    scalar nbZeta = curNearestPoint.second();

    // Info << "_: " << bI << ", " << segI << ", " << nbI << ", "
    // 	 << nbSegI << ", " << nbZeta << endl;
    
	// This set() takes the solver to lineContacts.C and executes the set() there
    lineContacts_[bI][segI][nbI].set
    (
        bI,
	segI,
	0,
	nbI,
	nbSegI,
	nbZeta,
	beam_.runTime().timeIndex()
    );
}


void Foam::beamContactModel::initializeLineContacts()
{
    // Info << "Inside initializeLineContacts() in beamContactModel.C" << endl;
    label nBeams = beam_.mesh().cellZones().size();

    // Update nearest points and contact angles
    
    label start = 0;
    
    for(label bI = 0; bI<nBeams; bI++)
    {
        for (label segI=0; segI<splines_[bI].nSegments(); segI++)
	{
            label globalSegIndex = start + segI;

            label cellID =
                 beam_.localCellIndex(globalSegIndex);

            if (cellID != -1)
            {
                for (label nbI=0; nbI<nBeams; nbI++)
                {
                    if (nbI != bI) // No self-contact
                    {
                        initializeLineContact(bI, segI, nbI);
                    }
                }
            }
	}

        start += splines_[bI].nSegments();
    }
}

void Foam::beamContactModel::updateLineContacts()
{
	// Info << "Inside updateLineContacts() beamContactModel.C" << endl;
    label nBeams = beam_.mesh().cellZones().size();

    // Update nearest points and contact angles
    label start = 0;
    for(label bI=0; bI<nBeams; bI++)
    {
	for (label segI=0; segI<splines_[bI].nSegments(); segI++)
	{
            label globalSegIndex = start + segI;

            label cellID =
                beam_.localCellIndex(globalSegIndex);

            if (cellID != -1)
            {
                for (label nbI=0; nbI<nBeams; nbI++)
                {
                    if (nbI != bI) // No self-contact
                    {
                        // initializeLineContact(bI, segI, nbI);
			
                        label oldNeiSegI =
                            lineContacts_[bI][segI][nbI].secondBeamSegment();

                        if (oldNeiSegI != -1)
                        {
                            updateLineContact(bI, segI, nbI);
                        }
                        else
                        {
                            // Info << "xxxxxxxxxxxxxxxx" << endl;
                            initializeLineContact(bI, segI, nbI);
                        }
                    }
                }
            }
            else
            {
                // Pout << "Global segment: " << globalSegIndex
                //      << " not present on this processor" << endl;
            }
	}

        start += splines_[bI].nSegments();
    }
}

void Foam::beamContactModel::updateLineContact
(
    const label bI,
    const label segI,
    const label nbI
)
{
    // Info << "updateLineContact: " <<  bI << ", " << segI << endl;

    const vectorField& neiPoints = splines_[nbI].points();
    const vectorField& neidRdS = splines_[nbI].dRdS();
    const vectorField& curCentres = splines_[bI].midPoints();

    label oldNeiSegI = lineContacts_[bI][segI][nbI].secondBeamSegment();

    if (oldNeiSegI == -1)
    {
        FatalErrorIn
		(
            "Foam::beamContactModel::updateLineContact(...) const"
		)
		<< "Trying to update line contact which is not initialized before."
		<< abort(FatalError);	    
    }

    label offset = 1;
    
    label fNeiSegI = oldNeiSegI-offset;
    
    if (fNeiSegI < 0)
    {
        fNeiSegI = 0;
    }

    label lNeiSegI = oldNeiSegI+offset;
    if (lNeiSegI > (splines_[nbI].nSegments()-1))
    {
        lNeiSegI = splines_[nbI].nSegments()-1;
    }
    
    labelScalar curNearestPoint(-1, 0);

    // Find nearest segment
    for
    (
        label nbSegI = fNeiSegI; nbSegI <= lNeiSegI; nbSegI++
    )
    {
        scalar testValue =
		(
			neidRdS[nbSegI]
			& (curCentres[segI] - neiPoints[nbSegI])
		)
		*(
			neidRdS[nbSegI + 1]
			& (curCentres[segI] - neiPoints[nbSegI + 1])
		);

		if (testValue <= 0)
		{
			curNearestPoint = 
				splines_[nbI].nearestPoint
				(
				nbSegI,
				curCentres[segI]
				);
			break;
		}
    }

    label neiSegI = curNearestPoint.first();

    if (neiSegI != -1)
    {
        scalar neiZeta = curNearestPoint.second();
	
        lineContacts_[bI][segI][nbI].set
	(
            bI,
	    segI,
	    0,
	    nbI,
	    neiSegI,
	    neiZeta,
	    beam_.runTime().timeIndex()
        );
    }
    else
    {
        // Try to initialize
        initializeLineContact(bI, segI, nbI);
    }
}


void Foam::beamContactModel::updateLineContactForces()
{
	// Info << "Inside updateLineContactForces() beamContactModel.C" << endl;
    label nBeams = beam_.mesh().cellZones().size();

    label start = 0;
    for(label bI=0; bI<nBeams; bI++)
    {
	for (label segI=0; segI<splines_[bI].nSegments(); segI++)
	{
            label globalSegIndex = start + segI;

            label cellID =
                beam_.localCellIndex(globalSegIndex);

            if (cellID != -1)
            {
                for (label nbI=0; nbI<nBeams; nbI++)
                {
                    if (nbI != bI) // No self-contact
                    {
                        // Pout << bI << ", " << segI << endl;

                        lineContacts_[bI][segI][nbI].updateForce
                        (
                            beam_,
                            splines_,
                            normalModel_(),
                            frictionModel_(),
			    lowerContactAngleLimit_,
                            upperContactAngleLimit_,
                            augmentedLagrangian_
                        );
                    }
                }
            }
	}
        
        start += splines_[bI].nSegments();
    }
}


void Foam::beamContactModel::updateNormalGapOffset()
{
    label nBeams = beam_.mesh().cellZones().size();

    // Calculate average penetration between two beams
    scalarFieldField avgPen(nBeams);
    scalarFieldField maxPen(nBeams);
    for(label bI=0; bI<nBeams; bI++)
    {
        avgPen.set(bI, new scalarField(nBeams, 0));
        maxPen.set(bI, new scalarField(nBeams, 0));

        for (label nbI=0; nbI<nBeams; nbI++)
        {
            if (nbI != bI) // No self-contact
            {
                scalar maxPenetration = 0;
                scalar avgPenetration = 0;
                scalar nContactSeg = 0;
                for (label segI=0; segI<splines_[bI].nSegments(); segI++)
                {
                    scalar curPenetration =
                       -lineContacts_[bI][segI][nbI].calcNormalGap
                        (
                            beam_,
                            splines_
                        );
                    
                    if (curPenetration > 0)
                    {
                        avgPenetration += curPenetration;
                        nContactSeg += 1;

                        if (curPenetration > maxPenetration)
                        {
                            maxPenetration = curPenetration;
                        }
                    }
                }

                avgPenetration /= nContactSeg + SMALL;
                
                avgPen[bI][nbI] = avgPenetration;
                maxPen[bI][nbI] = maxPenetration;
            }
        }

        // Info << bI << ", avgPen: " << avgPen[bI]
        //      << ", maxPen: " << maxPen[bI] << endl;
    }

    // // Update normal gap offset
    // for (label bI=0; bI<nBeams; bI++)
    // {
    //     for (label nbI=0; nbI<nBeams; nbI++)
    //     {
    //         if (nbI != bI) // No self-contact
    //         {
    //             scalar avg = (avgPen[bI][nbI] + avgPen[nbI][bI])/2;
                
    //             for (label segI=0; segI<splines_[bI].nSegments(); segI++)
    //             {
    //                 lineContacts_[bI][segI][nbI].normalGapOffset() += avg;
    //             }
    //         }
    //     }
    // }
}


void Foam::beamContactModel::initializePointContacts()
{
     // Info << "Initialize point contacts" << endl;
 
    for(label bI=0; bI<splines_.size(); bI++)
    {
        for (label nbI=bI; nbI<splines_.size(); nbI++)
        {
            if (nbI != bI) // No self-contact
            {
                for (label segI=0; segI<splines_[bI].nSegments(); segI++)
                {
                    // Find nearest segment
                    for
                    (
                        label neiSegI=0;
                        neiSegI<splines_[nbI].nSegments();
                        neiSegI++
                    )
                    {
                        // Check for point contact
                        Tuple2<scalar, scalar> zetas =
                            splines_[bI].checkPointContact
                            (
                                segI,
                                splines_[nbI],
                                neiSegI,
                                lowerContactAngleLimit_
                            );
                        scalar zetaC = zetas.first();
                        scalar neiZetaC = zetas.second();
			
			if 
			(
			    (zetaC > -1.5 && zetaC < 1.5) &&
			    (neiZetaC > -1.5 && neiZetaC < 1.5)
			)
			{
			    Info << "\nsegI " << segI << " zetaC "  << zetaC 
				 << "\nneiSegI " <<  neiSegI << " neiZetaC " << neiZetaC
				 << endl;
			}

                        if
                        (
			//- SB changed (18 Aug 2022): > -1 to >= -1
                            (zetaC > -1) && (zetaC <= 1)
                         && (neiZetaC > -1) && (neiZetaC <= 1)
                        )
                        {			     
                            // Check contact angle
                            vector tan0 =
                                splines_[bI].firstDerivative(segI, zetaC);
                            tan0 /= mag(tan0) + SMALL;

                            vector tan1 =
                                splines_[nbI].firstDerivative
                                (
                                    neiSegI,
                                    neiZetaC
                                );
                            tan1 /= mag(tan1) + SMALL;

                            scalar contactAngle =
                                ::acos(mag(tan0 & tan1))*180.0/M_PI;
                                
			    Info << "contact Angle " << tab << contactAngle << endl;
							
                            //if (contactAngle > lowerContactAngleLimit_)
                            //{
                                if (bI < nbI)
                                {
                                    pointContacts_.append
                                    (
                                        pointContact
                                        (
                                            bI,
                                            segI,
                                            zetaC,
                                            nbI,
                                            neiSegI,
                                            neiZetaC,
                                            beam_.runTime().timeIndex()
                                        )
                                    );
                                }
                                else
                                {
                                    pointContacts_.append
                                    (
                                        pointContact
                                        (
                                            nbI,
                                            neiSegI,
                                            neiZetaC,
                                            bI,
                                            segI,
                                            zetaC,
                                            beam_.runTime().timeIndex()
                                        )
                                    );
                                }
                            //}
                        }
                    }
                }
            }
        }
    }

    // Info << "Initialize point contacts: finished" << endl;
    
    // Calculate current point contact forces
    updatePointContactForces();
}

void Foam::beamContactModel::updatePointContacts()
{
    Info << "Inside updatePointContacts() beamContactModel.C" << endl;
    forAll(pointContacts_, pcI)
    {
        label bI = pointContacts_[pcI].firstBeam();
        label segI = pointContacts_[pcI].firstBeamSegment();
        // scalar zeta = pointContacts_[pcI].firstBeamZeta();

        label nbI = pointContacts_[pcI].secondBeam();
        label neiSegI = pointContacts_[pcI].secondBeamSegment();
        // scalar neiZeta = pointContacts_[pcI].secondBeamZeta();

        label fSegI = segI-2;
        if (fSegI < 0)
        {
            fSegI = 0;
        }
        label lSegI = segI+2;
        if (lSegI > (splines_[bI].nSegments()-1))
        {
            lSegI = splines_[bI].nSegments()-1;
        }
        
        label fNeiSegI = neiSegI-2;
        if (fNeiSegI < 0)
        {
            fNeiSegI = 0;
        }
        label lNeiSegI = neiSegI+2;
        if (lNeiSegI > (splines_[nbI].nSegments()-1))
        {
            lNeiSegI = splines_[nbI].nSegments()-1;
        }

        bool updated = false;
       
        for (label sI=fSegI; sI<=lSegI; sI++)
        {
            // Find nearest segment
            for (label nsI=fNeiSegI; nsI<=lNeiSegI; nsI++)
            {
                // Check for point contact
                Tuple2<scalar, scalar> zetas =
                    splines_[bI].checkPointContact
                    (
                        sI,
                        splines_[nbI],
                        nsI,
                        lowerContactAngleLimit_
                    );
                scalar zetaC = zetas.first();
                scalar neiZetaC = zetas.second();
		
		if 
		(
		    (zetaC > -1.5 && zetaC < 1.5) &&
		    (neiZetaC > -1.5 && neiZetaC < 1.5)
		)
		{
		    Info << "\nsegI " << sI << " updated zeta " << zetaC 
			 << "\nneiSegI " << nsI << " updated neiZeta " << neiZetaC
			 << endl;
		}

                if
                (
                    (zetaC > -1) && (zetaC <= 1)
                 && (neiZetaC > -1) && (neiZetaC <= 1)
                )
                {
                    // Check contact angle
                    vector tan0 =
                        splines_[bI].firstDerivative(sI, zetaC);
                    tan0 /= mag(tan0);

                    vector tan1 =
                        splines_[nbI].firstDerivative
                        (
                            nsI,
                            neiZetaC
                        );
                    tan1 /= mag(tan1);

                    scalar contactAngle =
                        ::acos(mag(tan0 & tan1))*180.0/M_PI;

                   // if (contactAngle > lowerContactAngleLimit_)
                   // {
                        if (bI < nbI)
                        {
                            pointContacts_[pcI].set
                            (
                                bI,
                                sI,
                                zetaC,
                                nbI,
                                nsI,
                                neiZetaC,
                                beam_.runTime().timeIndex()
                            );
                        }
                        else
                        {
                            pointContacts_[pcI].set
                            (
                                nbI,
                                nsI,
                                neiZetaC,
                                bI,
                                sI,
                                zetaC,
                                beam_.runTime().timeIndex()
                            );
                        }
                        
                        updated = true;
                        break;
                   // }
                }
            }
        }

        if (!updated)
        {
            Info << "Point contact " << pointContacts_[pcI]
                 << " is not updated." << endl;
        }
    }
    
    // Calculate current point contact forces
    updatePointContactForces();
}

void Foam::beamContactModel::addNewPointContacts()
{
    for (label bI=0; bI<splines_.size(); bI++)
    {
        for (label nbI=bI; nbI<splines_.size(); nbI++)
        {
            if (nbI != bI) // No self-contact
            {
                for (label segI=0; segI<splines_[bI].nSegments(); segI++)
                {
                    // Find nearest segment
                    for
                    (
                        label neiSegI=0;
                        neiSegI<splines_[nbI].nSegments();
                        neiSegI++
                    )
                    {
                        // Check for point contact
                        Tuple2<scalar, scalar> zetas =
                            splines_[bI].checkPointContact
                            (
                                segI,
                                splines_[nbI],
                                neiSegI,
                                lowerContactAngleLimit_
                            );
                        scalar zetaC = zetas.first();
                        scalar neiZetaC = zetas.second();
			
	    
                        if
                        (
                            (zetaC > -1) && (zetaC <= 1)
                         && (neiZetaC > -1) && (neiZetaC <= 1)
                        )
                        {
			    
                            // Check contact angle
                            vector tan0 =
                                splines_[bI].firstDerivative(segI, zetaC);
                            tan0 /= mag(tan0);

                            vector tan1 =
                                splines_[nbI].firstDerivative
                                (
                                    neiSegI,
                                    neiZetaC
                                );
                            tan1 /= mag(tan1);

                            scalar contactAngle =
                                ::acos(mag(tan0 & tan1))*180.0/M_PI;

                          //  if (contactAngle > lowerContactAngleLimit_)
                           // {
                                pointContact newPointContact;
                                
                                if (bI < nbI)
                                {
                                    newPointContact.set
                                    (
                                        bI,
                                        segI,
                                        zetaC,
                                        nbI,
                                        neiSegI,
                                        neiZetaC,
                                        beam_.runTime().timeIndex()
                                    );
                                }
                                else
                                {
                                    newPointContact.set
                                    (
                                        nbI,
                                        neiSegI,
                                        neiZetaC,
                                        bI,
                                        segI,
                                        zetaC,
                                        beam_.runTime().timeIndex()
                                    );
                                }
                                
                                // Check if this contact point
                                // already exists
                                bool exist = false;
                                forAll(pointContacts_, pcI)
                                {
                                    if (pointContacts_[pcI] == newPointContact)
                                    {
                                        exist = true;
                                        break;
                                    }
                                }

                                if (!exist)
                                {
                                     Info << "contactAngle: "
                                         << contactAngle << endl;
                                    
                                    newPointContact.updateForce
                                    (
                                        beam_,
                                        splines_,
                                        pointNormalModel_(),
                                        pointFrictionModel_(),
					lowerContactAngleLimit_,
					upperContactAngleLimit_,
                                        augmentedLagrangian_
                                    );

                                    if
                                    (
                                        newPointContact.delta()
                                      < 4*(beam_.R(bI)+beam_.R(nbI)/2)
                                    )
                                    {
                                        pointContacts_.append(newPointContact);
                                    }
                                     Info << "adding contact point: "
                                          << newPointContact << endl;
                                }
                          //  }
                        }
                    }
                }
            }
        }
    }
}

void Foam::beamContactModel::updatePointContactForces()
{
    // Calculate current point contact forces
    forAll(pointContacts_, pcI)
    {
        pointContacts_[pcI].updateForce
        (
            beam_,
            splines_,
            pointNormalModel_(),
            pointFrictionModel_(),
	    lowerContactAngleLimit_,
	    upperContactAngleLimit_,
            augmentedLagrangian_
        );
    }

     Info << "pointContacts_.size(): " << pointContacts_.size() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::beamContactModel::beamContactModel(const beamModel& beam)
:
    IOdictionary
    (
        IOobject
        (
            "beamContactProperties",
            beam.mesh().time().constant(),
            beam.mesh().time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    beam_(beam),
    splines_(),
    normalModel_
    (
        normalContactModel::New
        (
            //   word(lookup("normalContactModel")),
            //   *this
	    word(subDict("lineContact").lookup("normalContactModel")),
	    subDict("lineContact")
        )
    ),
    frictionModel_
    (
        frictionContactModel::New
        (
            //   word(lookup("frictionContactModel")),
            //   *this
	    word(subDict("lineContact").lookup("frictionContactModel")),
	    subDict("lineContact")
        )
    ),
    pointNormalModel_
    (
	normalContactModel::New
	(
	    word(subDict("pointContact").lookup("normalContactModel")),
	    subDict("pointContact")
	)
    ),
    pointFrictionModel_
    (
	frictionContactModel::New
	(
	    word(subDict("pointContact").lookup("frictionContactModel")),
	    subDict("pointContact")
	)
    ),
    lineContacts_(),
    conicalPulleyContacts_(),
    toroidalPulleyContacts_(),
    conicalPulleyLoads_(0),
    conicalPulleyContactGaps_(),
    activeConicalPulleyContact_
    (
        IOobject
        (
            "activeConicalPulleyContact",
            beam.runTime().timeName(),
            beam.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beam.mesh(),
        dimensionedScalar("zero", dimless, -1)
    ),
    conicalPulleyNormalContactForce_
    (
        IOobject
        (
            "conicalPulleyNormalContactForce",
            beam.runTime().timeName(),
            beam.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beam.mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    conicalPulleyAxialContactForce_
    (
        IOobject
        (
            "conicalPulleyAxialContactForce",
            beam.runTime().timeName(),
            beam.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beam.mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    conicalPulleyContactForceRatio_
    (
        IOobject
        (
            "conicalPulleyFctByFcn",
            beam.runTime().timeName(),
            beam.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beam.mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    toroidalPulleyLoads_(0),
    toroidalPulleyContactGaps_(),
    activeToroidalPulleyContact_
    (
        IOobject
        (
            "activeToroidalPulleyContact",
            beam.runTime().timeName(),
            beam.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beam.mesh(),
        dimensionedScalar("zero", dimless, -1)
    ),
    toroidalPulleyNormalContactForce_
    (
        IOobject
        (
            "toroidalPulleyNormalContactForce",
            beam.runTime().timeName(),
            beam.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beam.mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    toroidalPulleyContactForceRatio_
    (
        IOobject
        (
            "toroidalPulleyFctByFcn",
            beam.runTime().timeName(),
            beam.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        beam.mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    nearestPoints_(),
    contactForces_(),
    contactForceDerivatives_(),
    contactAngles_(),
    pointContacts_(),
    frictionalContactForces_(),
    frictionalContactMoments_(),
    contactGaps_(),
    contactOffsets_(),
    implicit_(lookupOrDefault<bool>("implicit", true)),
    augmentedLagrangian_(lookupOrDefault<bool>("augmentedLagrangian", false)),
	sendAllPointsAndTangents_
    (
        lookupOrDefault<bool>("sendAllPointsAndTangents", true)
    ),
    lowerContactAngleLimit_(lookupOrDefault<scalar>("lowerContactAngleLimit", 15)),
    upperContactAngleLimit_(lookupOrDefault<scalar>("upperContactAngleLimit", 30)),
    curTimeIndex_(-1),
    totalSplinesUpdateTime_(0)
{
    label nBeams = beam_.mesh().cellZones().size();
    
    Info << "lower contact angle " << tab << lowerContactAngleLimit_ << endl;
    Info << "upper contact angle " << tab << upperContactAngleLimit_ << endl;

    // Create splines
    {
        splines_.clear();
        splines_.setSize(nBeams);
        forAll(splines_, sI)
        {
            vectorField curPoints;
            vectorField curTangents;
            beam_.currentGlobalBeamPointsAndTangents
            (
                sI,
                curPoints,
                curTangents
            );

            splines_.set
            (
                sI,
                new HermiteSpline(curPoints, curTangents)
                // new HermiteSpline
                // (
                //     beam_.currentGlobalBeamPoints(sI),
                //     beam_.currentGlobalBeamTangents(sI)
                // )
            );
        }
    }

    // if (Pstream::myProcNo() == 0)
    // {
    //     Pout << "length[0]: "
    //          << splines_[0].length() << endl;
    //     Pout << "length[1]: "
    //          << splines_[1].length() << endl;
    // }

    // if (Pstream::myProcNo() == 1)
    // {
    //     Pout << "length[0]: "
    //          << splines_[0].length() << endl;
    //     Pout << "length[1]: "
    //          << splines_[1].length() << endl;
    // }

    // Initialiye line contacts
    lineContacts_.setSize(nBeams);
    forAll(lineContacts_, bI)
    {
        lineContacts_.set
        (
            bI,
            new lineContactListList
            (
                splines_[bI].nSegments()
            )
        );

        forAll(lineContacts_[bI], segI)
        {
	    lineContacts_[bI][segI] =
	        lineContactList
	        (
		    nBeams,
		    lineContact()
                );
        }
    }

    initializeLineContacts();

    // // Print line contacts
    // for(label bI=0; bI<nBeams; bI++)
    // {
    //     for (label segI=0; segI<splines_[bI].nSegments(); segI++)
    //     {
    //         for (label nbI=0; nbI<nBeams; nbI++)
    //         {
    //             if (nbI != bI) // No self-contact
    //             {
    //                 if (Pstream::myProcNo() == 1)
    //                 {
    //                     label nbSegI =
    //                         lineContacts_[bI][segI][nbI]
    //                        .secondBeamSegment();

    //                     Pout << bI << ", " << segI << ", "
    //                          << nbI << ", " << nbSegI << endl;
    //                 }
    //     	}
    //         }
    //     }
    // }

    // Pout << "print initial contact" << endl;
    // sleep(5);

    // Initialiye conical pulley contacts
    label nPulleys = beam_.conicalPulleys().size();
    if (nPulleys)
    {
        conicalPulleyContacts_.setSize(nBeams);

        conicalPulleyContactGaps_.setSize(nBeams);

        forAll(conicalPulleyContacts_, bI)
        {
            conicalPulleyContacts_.set
            (
                bI,
                new conicalPulleyContactListList
                (
                    splines_[bI].nSegments()
                )
            );

            conicalPulleyContactGaps_.set
            (
                bI,
                new scalarPairList(nPulleys, scalarPair(GREAT, GREAT))
                // new scalarField(nPulleys, GREAT)
            );

            forAll(conicalPulleyContacts_[bI], segI)
            {
                conicalPulleyContacts_[bI][segI] =
                    conicalPulleyContactList
                    (
                        nPulleys,
                        conicalPulleyContact(conicalPulleyContacts_)
                    );
            }
        }

        updateConicalPulleyContacts();

        conicalPulleyLoads_.setSize(nPulleys, vector::zero);
    }


    // Initialiye toroidal pulley contacts
    nPulleys = beam_.toroidalPulleys().size();
    if (nPulleys)
    {
        toroidalPulleyContacts_.setSize(nBeams);

        toroidalPulleyContactGaps_.setSize(nBeams);

        forAll(toroidalPulleyContacts_, bI)
        {
            toroidalPulleyContacts_.set
            (
                bI,
                new toroidalPulleyContactListList
                (
                    splines_[bI].nSegments()
                )
            );

            toroidalPulleyContactGaps_.set
            (
                bI,
                new scalarPairList(nPulleys, scalarPair(GREAT, GREAT))
                // new scalarField(nPulleys, GREAT)
            );

            forAll(toroidalPulleyContacts_[bI], segI)
            {
                toroidalPulleyContacts_[bI][segI] =
                    toroidalPulleyContactList
                    (
                        nPulleys,
                        toroidalPulleyContact()
                    );
            }
        }

        updateToroidalPulleyContacts();

        toroidalPulleyLoads_.setSize(nPulleys, vector::zero);
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
                        splines_[bI].nSegments(),
                        labelScalar(0, 0)
                    );
            }
        }
    }

    // Initialize contact angles
    contactAngles_.setSize(nBeams);
    forAll(contactAngles_, bI)
    {
        contactAngles_.set
        (
            bI,
            new scalarListList
            (
                nBeams
            )
        );

        forAll(contactAngles_[bI], nbI)
        {
            if (nbI != bI)
            {
                contactAngles_[bI][nbI] =
                    List<scalar>
                    (
                        splines_[bI].nSegments(),
                        scalar(0)
                    );
            }
        }
    }

    // Initialize contact forces
    contactForces_.setSize(nBeams);
    contactForceDerivatives_.setSize(nBeams);
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
        contactForceDerivatives_.set
        (
            bI,
            new tensorListList
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
                        splines_[bI].nSegments(),
                        vector::zero
                    );
                contactForceDerivatives_[bI][nbI] =
                    tensorList
                    (
                        splines_[bI].nSegments(),
                        tensor::zero
                    );
                frictionalContactForces_[bI][nbI] =
                    vectorList
                    (
                        splines_[bI].nSegments(),
                        vector::zero
                    );
                frictionalContactMoments_[bI][nbI] =
                    vectorList
                    (
                        splines_[bI].nSegments(),
                        vector::zero
                    );
            }
        }
    }

    // Initialize contact gap
    contactGaps_.setSize(nBeams);
    forAll(contactGaps_, bI)
    {
        contactGaps_.set
        (
            bI,
            new scalarListList
            (
                nBeams
            )
        );
        forAll(contactGaps_[bI], nbI)
        {
            if (nbI != bI)
            {
                contactGaps_[bI][nbI] =
                    scalarList
                    (
                        splines_[bI].nSegments(),
                        GREAT
                    );
            }
        }
    }

    // Initialize contact offset
    contactOffsets_.setSize(nBeams);
    forAll(contactOffsets_, bI)
    {
        contactOffsets_.set
        (
            bI,
            new scalarListList
            (
                nBeams
            )
        );
        forAll(contactOffsets_[bI], nbI)
        {
            if (nbI != bI)
            {
                contactOffsets_[bI][nbI] =
                    scalarList
                    (
                        splines_[bI].nSegments(),
                        0.0
                    );
            }
        }
    }

    // Initialize contact distance
    contactDistances_.setSize(nBeams);
    forAll(contactDistances_, bI)
    {
        contactDistances_.set
        (
            bI,
            new scalarListList
            (
                nBeams
            )
        );
        forAll(contactDistances_[bI], nbI)
        {
            if (nbI != bI)
            {
                contactDistances_[bI][nbI] =
                    scalarList
                    (
                        splines_[bI].nSegments(),
                        GREAT
                    );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::beamContactModel::~beamContactModel()
{
    // deleteDemandDrivenData(normalModelPtr_);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::beamContactModel::update()
{
    Info << "Inside update() function of beamContactModel.C file" << endl;
    bool debug
    (
        lookupOrDefault<bool>
        (
            "debug",
            false
        )
    );

    label nBeams = beam_.mesh().cellZones().size();
    
    // Info << "nBeams: \n " << nBeams << endl;

    // Update spline geometry
    if (beam_.iOuterCorr() == 0)
    {
        // All points and tangents have been updated at the end
        // of the previous time step
    }
    else if (sendAllPointsAndTangents_)
    {
        scalar tStart = beam_.runTime().elapsedCpuTime();
        
        // Info << "tStart: \n" << tStart << endl;

        // Update splines
        forAll(splines_, sI)
        {
            vectorField curPoints;
            vectorField curTangents;
	    // Info << "I am inside update splines \n" << endl;
            // Transfer all points
            beam_.currentGlobalBeamPointsAndTangents
            (
                sI,
                curPoints,
                curTangents
                );

            splines_[sI].movePoints(curPoints, curTangents);
        }

        scalar tEnd = beam_.runTime().elapsedCpuTime();

        totalSplinesUpdateTime_ += tEnd - tStart;
        
    }
    else
    {
        // Update processor points and tangents and
        // only active (in-contact) off-processor points and tangents

        // First identify spline points beloniging to this processor,
        // which have to be send to the neighbour processors (active points)
        
        
	// Info << "I am in else loop: \n" << endl;
		
        scalar tStart = beam_.runTime().elapsedCpuTime();
        
        PtrList<List<labelHashSet> > splPtsFromProc(splines_.size());
        for(label bI=0; bI<nBeams; bI++)
        {
            splPtsFromProc.set
            (
                bI,
                new List<labelHashSet>(Pstream::nProcs(), labelHashSet())
            );
        }

        // label start = 0;
        for(label bI=0; bI<nBeams; bI++)
        {
            for (label segI=0; segI<splines_[bI].nSegments(); segI++)
            {
                label globalSegIndex =
                    beam_.startCells()[bI] + segI;

                label locCellID = beam_.localCellIndex(globalSegIndex);

                if (locCellID != -1)
                {
                    for (label nbI=0; nbI<nBeams; nbI++)
                    {
                        if (nbI != bI) // No self-contact
                        {
                            labelList glNeiCellIDs(2, -1);

                            glNeiCellIDs[0] =
                                beam_.startCells()[nbI]
                              + lineContacts_[bI][segI][nbI].
                                secondBeamLowerSegment();

                            glNeiCellIDs[1] =
                                beam_.startCells()[nbI]
                              + lineContacts_[bI][segI][nbI].
                                secondBeamUpperSegment();

                            forAll(glNeiCellIDs, ncI)
                            {
                                label locNeiCellId =
                                    beam_.localCellIndex(glNeiCellIDs[ncI]);

                                if (locNeiCellId == -1)
                                {
                                    // Pout << glNeiCellIDs << endl;
                            
                                    label glNeiSegID =
                                        beam_.whichSegment(glNeiCellIDs[ncI]);

                                    labelPair procLocCellID =
                                        beam_.procLocalCellIndex
                                        (
                                            glNeiCellIDs[ncI]
                                        );

                                    label neiProcID = procLocCellID.first();
                                    // label locCellID = procLocCellID.second();

                                    label firstGlPoint = glNeiSegID;
                                    label secondGlPoint = glNeiSegID+1;
                                    
                                    if
                                    (
                                       !splPtsFromProc[nbI][neiProcID].found
                                        (
                                            firstGlPoint
                                        )
                                    )
                                    {
                                        splPtsFromProc[nbI][neiProcID].
                                            insert(firstGlPoint);
                                    }

                                    if
                                    (
                                       !splPtsFromProc[nbI][neiProcID].found
                                        (
                                            secondGlPoint
                                        )
                                    )
                                    {
                                        splPtsFromProc[nbI][neiProcID].
                                            insert(secondGlPoint);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Send global point data to corresponding processors
        PtrList<labelListList> splinePointsSendToProc_(nBeams);
        for (label bI=0; bI<nBeams; bI++)
        {
            splinePointsSendToProc_.set
            (
                bI,
                new labelListList(Pstream::nProcs())
            );

            for (label procI=0; procI<Pstream::nProcs(); procI++)
            {
                if (procI != Pstream::myProcNo())
                {
                    labelList curPts =
                        splPtsFromProc[bI][procI].toc();

                    {
                        OPstream toProc
                        (
                            Pstream::blocking,
                            procI,
                            sizeof(label)
                        );
                        
                        toProc << curPts.size();
                    }
                    
                    {   
                        // Parallel data exchange
                        OPstream::write
                        (
                            Pstream::blocking,
                            procI,
                            reinterpret_cast<const char*>
                            (
                                curPts.begin()
                            ),
                            curPts.byteSize()
                        );
                    }                    
                }
            }

            for (label procI=0; procI<Pstream::nProcs(); procI++)
            {
                if (procI != Pstream::myProcNo())
                {
                    label numOfPoints = 0;
                    {
                        IPstream fromProc
                        (
                            Pstream::blocking,
                            procI,
                            sizeof(label)
                        );
                        
                        fromProc >> numOfPoints;
                    }

                    {
                        labelList& curPoints =
                            splinePointsSendToProc_[bI][procI];
                        curPoints.setSize(numOfPoints);

                        // Parallel data exchange
                        IPstream::read
                        (
                            Pstream::blocking,
                            procI,
                            reinterpret_cast<char*>
                            (
                                curPoints.begin()
                            ),
                            curPoints.byteSize()
                        );

                        // Transform global to local indices
                        forAll(curPoints, pI)
                        {
                            label tmpGlPtIndex = curPoints[pI];
                            
                            curPoints[pI] =
                                beam_.globalToLocalBeamPointsAddressing()
                                [bI][curPoints[pI]];

                            if (curPoints[pI] == -1)
                            {
                                FatalErrorIn
                                (
                                    "Foam::beamContactModel::update() const"
                                )
                                    << "Global point index: "
                                    << tmpGlPtIndex << " of beam: " << bI
                                    << " does not correspond to any of the"
                                    << " local point indices"
                                    << abort(FatalError);
                            }
                        }
                    }                    
                }
            }
        }
    
        forAll(splines_, sI)
        {
            vectorField curPoints = splines_[sI].points();
            vectorField curTangents = splines_[sI].tangents();

            beam_.currentGlobalBeamPointsAndTangents
            (
                sI,
                splinePointsSendToProc_[sI],
                curPoints,
                curTangents
            );

            splines_[sI].movePoints(curPoints, curTangents);
        }

        scalar tEnd = beam_.runTime().elapsedCpuTime();

        totalSplinesUpdateTime_ += tEnd - tStart;
    }

    if (debug)
    {
        Info << "Total spline update time: \n "
             << totalSplinesUpdateTime_ << endl;
    }
    
    if (debug)
    {
        Info << "Update conical pulley contact forces." << endl;
    }

    // Update conical pulleys contact
    updateConicalPulleyContacts();
    updateConicalPulleyContactForces();

    if (debug)
    {
        Info << "Update toroidal pulley contact forces." << endl;
    }

    // Update toroidal pulleys contact
    updateToroidalPulleyContacts();
    updateToroidalPulleyContactForces();

    if (debug)
    {
        Info << "Contact points calculation." << endl;
    }

    updateLineContacts();
    
    // if (curTimeIndex_ < beam_.runTime().timeIndex())
    // {
    //     updateNormalGapOffset();
    // }

    updateLineContactForces();

    // Update nearest points and contact angles
    if (false) // SB: changed this to true to check whether point contact is working, 04/05/2022
    {
    for(label bI=0; bI<nearestPoints_.size(); bI++)
    {
        labelScalarListList& curNearestPoints = nearestPoints_[bI];
        scalarListList& curContactAngles = contactAngles_[bI];
        const vectorField& curCentres = splines_[bI].midPoints();

        for (label nbI=0; nbI<nearestPoints_.size(); nbI++)
        {
            if (nbI != bI) // No self-contact
            {
                const vectorField& neiPoints = splines_[nbI].points();
                const vectorField& neidRdS = splines_[nbI].dRdS();

                DynamicList<labelScalarLabelScalar> curPointContactPoints;
                
                for (label segI=0; segI<curNearestPoints[nbI].size(); segI++)
                {
                    // Find nearest segment
                    for
                    (
                        label nbSegI=0;
                        nbSegI<splines_[nbI].nSegments();
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
                            curNearestPoints[nbI][segI] = 
                            // nearestPoints_[bI][nbI][segI] = 
                                splines_[nbI].nearestPoint
                                (
                                    nbSegI,
                                    curCentres[segI]
                                );
                            break;
                        }
                    }

                    // Calculating potential contact angle
                    vector tan0 = splines_[bI].firstDerivative(segI, 0);
                    tan0 /= mag(tan0)+SMALL;
                    
                    label nbSegI = nearestPoints_[bI][nbI][segI].first();
                    label nbZeta = nearestPoints_[bI][nbI][segI].second();
                    vector tan1 = splines_[nbI].firstDerivative(nbSegI, nbZeta);
                    tan1 /= mag(tan1)+SMALL;

                    curContactAngles[nbI][segI] =
                       ::acos(mag(tan0 & tan1))*180.0/M_PI;
                }

                if (debug)
                {
                    Info << "Contact angle (" << bI << ", " << nbI << "), avg: "
                         << average(curContactAngles[nbI])
                         << ", max: " << max(curContactAngles[nbI])
                         << ", min: " << min(curContactAngles[nbI]) << endl;
                }
            }
        }
    }
    }

    if (debug)
    {
        Info << "Calculating point contacts" << endl;
    }

    // // Calculating point contact
     bool reusePrevPointContacts
     (
         lookupOrDefault<bool>
         (
             "reusePrevPointContacts",
             false
         )
     );

    // Here the point contact initialization is done. 
    // The else loop is executed here
     if (reusePrevPointContacts && pointContacts_.size())
     {
         updatePointContacts();
        
         if (curTimeIndex_ < beam_.runTime().timeIndex())
         {
	    addNewPointContacts();
         }
     }
     else
     {
         pointContacts_.clear();
	 
	 Info << "Fresh point contact points calculated: beamContactModel.C" << endl;
         initializePointContacts();

         // Info << "Num of point contacts: "
         //      << pointContacts_.size() << endl;
     }

    if (debug)
    {
        Info << "Contact points calculation finished. "
           << "Calculating contact forces" << endl;
    }
    
    // Calculate contact force
    // scalar maxFc = maxContactForce();
    if (false) // SB: changed this to true to check whether point contact is working, 04/05/2022
    {
    scalar maxResidual = 0;
    for(label bI=0; bI<splines_.size(); bI++)
    {
        const vectorField& C = splines_[bI].midPoints();

        forAll(nearestPoints_[bI], nbI)
        {
            const labelScalarList& curNearestPoints = nearestPoints_[bI][nbI];
            const scalarList& curContactAngles = contactAngles_[bI][nbI];

            vectorList& curContactForces = contactForces_[bI][nbI];
            tensorList& curContactForceDerivatives =
                contactForceDerivatives_[bI][nbI];

            scalarList& curContactGaps = contactGaps_[bI][nbI];
            scalarList& curContactOffsets = contactOffsets_[bI][nbI];
            scalarList& curContactDistances = contactDistances_[bI][nbI];

            if (bI != nbI) // No self-contact
            {
                forAll(curNearestPoints, segI)
                {
                    vector curPoint =
                        splines_[nbI].position
                        (
                            curNearestPoints[segI].first(),
                            curNearestPoints[segI].second()
                        );

                    scalar dist = mag(C[segI] - curPoint);

                    curContactDistances[segI] = dist;
                    
                    vector n = C[segI] - curPoint;
                    n /= mag(n) + SMALL;

                    scalar gap = dist - beam_.R(bI) - beam_.R(nbI);

                    scalar curGapResidual =
                        2*mag(gap - curContactGaps[segI])
                       /(beam_.R(bI) + beam_.R(nbI));

                    curContactGaps[segI] = gap;

                    // Not line contact
                    if (curContactAngles[segI] > lowerContactAngleLimit_)
                    {
                        curContactForces[segI] = vector::zero;
                        curContactForceDerivatives[segI] = tensor::zero;
                    }
                    else
                    {
                        curContactForces[segI] =
                            normalModel_().contactForce
                            (
                                curContactGaps[segI]
                              - curContactOffsets[segI]
                            )*n;

                        curContactForceDerivatives[segI] =
                            normalModel_().contactForceDerivative
                            (
                                curContactGaps[segI]
                              - curContactOffsets[segI]
                            )*sqr(n);
                    }

                    if (curGapResidual > maxResidual)
                    {
                        maxResidual = curGapResidual;
                    }
                }
            }
        }
    }
    }

    if (debug)
    {
        // Nearest points
        for(label bI=0; bI<nearestPoints_.size(); bI++)
        {
            labelScalarListList& curNearestPoints = nearestPoints_[bI];

            for (label nbI=0; nbI<nearestPoints_.size(); nbI++)
            {
                if (nbI != bI) // No self-contact
                {
                    scalarField param(curNearestPoints[nbI].size(), 0);
                    for
                    (
                        label segI=0;
                        segI<curNearestPoints[nbI].size();
                        segI++
                    )
                    {
                        param[segI] = curNearestPoints[nbI][segI].second();
                    }

                    Info << "Nearest point segment parameters" << endl;
                    Info << "bI = " << bI << ", nbI = " << nbI << endl;
                    Info << "max: " << max(param) << ", avg: "
                         << average(param) << ", min: " << min(param) << endl;
                }
            }
        } 
    }

    if (curTimeIndex_ < beam_.runTime().timeIndex())
    {
        curTimeIndex_ = beam_.runTime().timeIndex();
    }
    
    return 0;
}


bool Foam::beamContactModel::active() const
{
    return (beam_.mesh().cellZones().size() > 1);
}


bool Foam::beamContactModel::finalUpdate()
{
    // Update splines
    forAll(splines_, sI)
    {
        vectorField curPoints;
        vectorField curTangents;

        // Transfer all points
        beam_.currentGlobalBeamPointsAndTangents
        (
            sI,
            curPoints,
            curTangents
        );

        splines_[sI].movePoints(curPoints, curTangents);
    }

    label nBeams = beam_.mesh().cellZones().size();
    label nPulleys = beam_.conicalPulleys().size();

    if (nPulleys)
    {
        // Update angular velocity of pulleys
        {
            label start = 0;
            scalarField omega(nPulleys, 0);
            labelList nContactSegments(nPulleys, 0);
            const vectorField& DW = beam_.solutionDW();
            // const vectorField& oldW = beam_.solutionW().oldTime();
            for(label bI=0; bI<nBeams; bI++)
            {
                for (label segI=0; segI<splines_[bI].nSegments(); segI++)
                {
                    label globalSegI = start + segI;
        
                    for (label pulleyI=0; pulleyI<nPulleys; pulleyI++)
                    {                        
                        scalarField magOldFcn(2, 0);
                        
                        magOldFcn[0] =
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .oldNormalContactForce()[0]
                            );
                    
                        magOldFcn[1] =
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .oldNormalContactForce()[1]
                            );
                        
                        // scalar maxOldFcn = max(magOldFcn);
                        // scalar minOldFcn = min(magOldFcn);
                        // scalar DOldFcn = maxOldFcn - minOldFcn;

                        scalarField magFcn(2, 0);
                    
                        magFcn[0] =
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .normalContactForce()[0]
                            );
                    
                        magFcn[1] =
                            mag
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .normalContactForce()[1]
                            );

                        // scalar maxFcn = max(magFcn);
                        // scalar minFcn = min(magFcn);
                        // scalar DFcn = maxFcn - minFcn;

                        scalar smallFcn =
                            conicalPulleyContacts_[bI][segI][pulleyI].smallFcn_;

                        if
                        (
                            (
                                (magOldFcn[0] > smallFcn)
                             || (magOldFcn[1] > smallFcn)
                            )
                         && (
                                (magFcn[0] > smallFcn)
                             || (magFcn[1] > smallFcn)
                            )
                        )
                        {
                            const vector& C = splines()[bI].midPoints()[segI];
                            vector T = splines()[bI].dRdS()[segI];
                            T /= mag(T) + SMALL;

                            vector DWa =
                                T*(T&DW[globalSegI]);
                                // T*(T&(W[globalSegI] - oldW[globalSegI]));

                            vector r =
                                C - beam_.conicalPulleys()[pulleyI].origin();
                            scalar sgn =
                                sign
                                (
                                    (r ^ DWa) &
                                    beam_.conicalPulleys()[pulleyI].axis()
                                );

                            // scalar phi =
                            // (
                            //     conicalPulleyContacts_[bI][segI][pulleyI]
                            //    .pulleyContactCoordinates()[0].y()
                            //   + conicalPulleyContacts_[bI][segI][pulleyI]
                            //    .pulleyContactCoordinates()[1].y()
                            // )/2;
                            
                            // scalar oldPhi =
                            // (
                            //     conicalPulleyContacts_[bI][segI][pulleyI]
                            //    .oldPulleyContactCoordinates()[0].y()
                            //   + conicalPulleyContacts_[bI][segI][pulleyI]
                            //    .oldPulleyContactCoordinates()[1].y()
                            // )/2;

                            // scalar DPhi = phi - oldPhi;

                            scalar avgContactRadius =
                            (
                                conicalPulleyContacts_[bI][segI][pulleyI]
                               .pulleyContactCoordinates()[0].x()
                              + conicalPulleyContacts_[bI][segI][pulleyI]
                               .pulleyContactCoordinates()[1].x()
                            )/2;

                            omega[pulleyI] +=
                                sgn*mag(DWa)/
                                beam_.conicalPulleys()[pulleyI].minRadius();
                                // avgContactRadius;
                                
                            
                            nContactSegments[pulleyI] += 1;
                        }

                        // 2100
                        // omega[pulleyI] = -2400/beam_.conicalPulleys()[pulleyI].minRadius();
                    }
                }

                start++;
            }

            forAll(omega, pulleyI)
            {
                // scalar Omega = omega[pulleyI];
                scalar Omega =
                    (omega[pulleyI]/(nContactSegments[pulleyI] + SMALL))
                   /beam_.runTime().deltaT().value();

                const_cast<conicalPulley&>(beam_.conicalPulleys()[pulleyI])
               .setAngularVelocity(Omega); //-2;

                Info << "Pulley: " << pulleyI << ", angular velocity: "
                     << beam_.conicalPulleys()[pulleyI].angularVelocity()
                     << " (" << Omega << ")" << endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
