/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

Class
    plasticityStressResultantReturn

\*---------------------------------------------------------------------------*/

#include "plasticityStressResultantReturn.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<plasticityStressResultantReturn> plasticityStressResultantReturn::New
(
    const word& name,
    beamModel& beamModel
)
{
    Info<< "\tPlasticity stress resultant return method: " << name << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(name);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "plasticityStressResultantReturn::New(\n"
            "    const word& name\n"
            "    beamModel& beamModel\n"
            ")",
            "beamModelCoeffs"
        )
            << "Unknown plasticityStressResultantReturn type "
            << name << endl << endl
            << "Valid  plasticityStressResultantReturns methods are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<plasticityStressResultantReturn>
    (
        cstrIter()
        (
            name,
            beamModel
        )
    );
}
    

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
