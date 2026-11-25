/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "lookupBeamModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const beamModel& lookupBeamModel(const objectRegistry& obReg)
{
    return lookupBeamModel(obReg, obReg.name());
}


const beamModel& lookupBeamModel
(
    const objectRegistry& obReg,
    const word& regionName
)
{
    // Beam models are registered on Time under a name that depends on the
    // order in which they were created ("beamProperties" for the first, see
    // the beamModel constructor), so match them by their mesh region
    wordList regions;

    for (const objectRegistry* reg : {&obReg, &obReg.parent()})
    {
        const HashTable<const beamModel*> beams =
            reg->lookupClass<beamModel>();

        forAllConstIters(beams, iter)
        {
            if (iter.val()->mesh().name() == regionName)
            {
                return *iter.val();
            }

            regions.append(iter.val()->mesh().name());
        }
    }

    FatalErrorInFunction
        << "Could not find the beamModel of region " << regionName << nl
        << "beamModels found, by region: " << regions
        << abort(FatalError);

    // Keep the compiler happy
    return obReg.lookupObject<beamModel>("none");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
