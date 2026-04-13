/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::functionObjects::beamContactData

Description
    FunctionObject reports BEAM contact data

SourceFiles
    beamContactData.C

\*---------------------------------------------------------------------------*/

#ifndef beamContactData_H
#define beamContactData_H

#include "fvMeshFunctionObject.H"
#include "OFstream.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class beamContactData Declaration
\*---------------------------------------------------------------------------*/

class beamContactData
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Contact history file pointer
        autoPtr<OFstream> historyFilePtr_;

    // Private Member Functions

        //- Write data
        bool writeData();

public:

    //- Runtime type information
    TypeName("beamContactData");


    // Constructors

        //- Construct from Time and dictionary
        beamContactData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        beamContactData(const beamContactData&) = delete;

        //- No copy assignment
        void operator=(const beamContactData&) = delete;


    //- Destructor
    virtual ~beamContactData() = default;


    // Member Functions

        //- Read the beamContactData data
        virtual bool read(const dictionary& dict);

        //- Execute the beamContactData
        virtual bool execute();

        //- Execute at the final time-loop
        virtual bool end();

        //- Write callback
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
