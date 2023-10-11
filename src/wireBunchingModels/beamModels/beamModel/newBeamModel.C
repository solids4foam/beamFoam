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

\*---------------------------------------------------------------------------*/

#include "beamModel.H"
//#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::beamModel> Foam::beamModel::New
(
    Time& runTime,
    const word& region
)
{
    word beamModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the flow is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary beamProperties
        (
            IOobject
            (
                "beamProperties",
                bool(region == dynamicFvMesh::defaultRegion)
              ? fileName(runTime.caseConstant())
              : fileName(runTime.caseConstant()/region),
                // runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        beamProperties.lookup("beamModel")
            >> beamModelTypeName;
    }

    Info<< "Selecting beamModel " << beamModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(beamModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "beamModel::New(Time&, const word&)"
        )   << "Unknown beamModel type " << beamModelTypeName
            << endl << endl
            << "Valid beamModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<beamModel>(cstrIter()(runTime, region));
}


// ************************************************************************* //
