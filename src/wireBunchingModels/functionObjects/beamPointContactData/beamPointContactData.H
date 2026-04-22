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
    Foam::functionObjects::beamPointContactData

Description
    FunctionObject reports normal point-contact data for each point contact.

SourceFiles
    beamPointContactData.C

\*---------------------------------------------------------------------------*/

#ifndef beamPointContactData_H
#define beamPointContactData_H

#include "fvMeshFunctionObject.H"
#include "OFstream.H"
#include "dictionary.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class beamPointContactData Declaration
\*---------------------------------------------------------------------------*/

class beamPointContactData
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

        //- Contact history directory
        fileName historyDir_;

        //- Contact history file pointers
        PtrList<OFstream> historyFilePtrs_;

    // Private Member Functions

        //- Create the history directory
        void makeHistoryDir();

        //- Create a point-contact history file
        void makeHistoryFile(const label pointContactI);

        //- Write data
        bool writeData();

public:

    //- Runtime type information
    TypeName("beamPointContactData");


    // Constructors

        //- Construct from Time and dictionary
        beamPointContactData
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        beamPointContactData(const beamPointContactData&) = delete;

        //- No copy assignment
        void operator=(const beamPointContactData&) = delete;


    //- Destructor
    virtual ~beamPointContactData() = default;


    // Member Functions

        //- Read the beamPointContactData data
        virtual bool read(const dictionary& dict);

        //- Execute the beamPointContactData
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
