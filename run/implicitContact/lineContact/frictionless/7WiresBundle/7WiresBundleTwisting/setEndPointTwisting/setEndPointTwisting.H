/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    setEndPointTwisting

Description
    FunctionObject set displacement boundary condition 
    at the beam end point

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

SourceFiles
    setEndPointTwisting.C

\*---------------------------------------------------------------------------*/

#ifndef setEndPointTwisting_H
#define setEndPointTwisting_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setEndPointTwisting Declaration
\*---------------------------------------------------------------------------*/

class setEndPointTwisting
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Region name
        word regionName_;

        //- Start patch index
        label startPatchIndex_;
  
        //- End patch index
        label endPatchIndex_;

        //- Beam end point patch (circle) radius
        scalar radius_;

        //- Twisting angle
        scalar thetaEnd_;
  
    // Private Member Functions

        //- Disallow default bitwise copy construct
        setEndPointTwisting
        (
            const setEndPointTwisting&
        );

        //- Disallow default bitwise assignment
        void operator=(const setEndPointTwisting&);

        //- Set boundary condition
        bool setBC();

public:

    //- Runtime type information
    TypeName("setEndPointTwisting");


    // Constructors

        //- Construct from components
        setEndPointTwisting
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
#if FOAMEXTEND > 40
        virtual bool execute(const bool forceWrite);
#else
        virtual bool execute();
#endif

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);
    
        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
