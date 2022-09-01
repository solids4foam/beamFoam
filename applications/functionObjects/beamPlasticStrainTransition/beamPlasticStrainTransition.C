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

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "beamPlasticStrainTransition.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(beamPlasticStrainTransition, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        beamPlasticStrainTransition,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::beamPlasticStrainTransition::updatePlasticStrain()
{
    Info << "Updating plastic strain" << endl;

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    surfaceVectorField& GammaP =
        const_cast<surfaceVectorField&>
        (
            mesh.lookupObject<surfaceVectorField>("GammaP")
        );

    surfaceVectorField& KP =
        const_cast<surfaceVectorField&>
        (
            mesh.lookupObject<surfaceVectorField>("KP")
        );

    scalar tStart = time_.startTime().value();    
    scalar tEnd = time_.endTime().value();
    scalar t = time_.value() + time_.deltaT().value();
    Info << tStart << ", " << tEnd << ", " << t << endl;

    GammaP = (*GammaPptr_)*(t - tStart)/(tEnd - tStart);

    KP = (*KPptr_)*(t - tStart)/(tEnd - tStart);

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::beamPlasticStrainTransition::beamPlasticStrainTransition
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    GammaPptr_(NULL),
    KPptr_(NULL)
    // historyFilePtr_(NULL),
    // patchName_(dict.lookup("patchName")),
    // patchIndex_(-1)
{
    Info << "Creating functio object " << name_ << "\n" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    GammaPptr_ =
        new surfaceVectorField
        (
            IOobject
            (
                "GammaP",
                time_.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

    KPptr_ =
        new surfaceVectorField
        (
            IOobject
            (
                "KP",
                time_.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

    // updatePlasticStrain();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::beamPlasticStrainTransition::start()
{
    Info << "start" << endl;
    return updatePlasticStrain();
    // return true;
}


#if FOAMEXTEND > 40
bool Foam::beamPlasticStrainTransition::execute(const bool forceWrite)
#else
bool Foam::beamPlasticStrainTransition::execute()
#endif
{
    Info << "execute" << endl;
    return updatePlasticStrain();
    // return true;
}


bool Foam::beamPlasticStrainTransition::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}

// ************************************************************************* //
