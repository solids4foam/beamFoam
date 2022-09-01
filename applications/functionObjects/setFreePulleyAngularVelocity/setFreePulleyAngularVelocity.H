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

Class
    setFreePulleyAngularVelocity

Description
    FunctionObject reports average traction at specified patch

Author
    Zeljko Tukovic, FSB.  All rights reserved

SourceFiles
    setFreePulleyAngularVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef setFreePulleyAngularVelocity_H
#define setFreePulleyAngularVelocity_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setFreePulleyAngularVelocity Declaration
\*---------------------------------------------------------------------------*/

class setFreePulleyAngularVelocity
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

        //- Pulley index
        labelList pulleyIndices_;
    
        //- Time when frictionless is changed
        scalar changeFrictionTime_;
    
    // Private Member Functions

        //- Set angular velocity for specified pulleys
        void setAngularVelocity();
    
        //- Disallow default bitwise copy construct
        setFreePulleyAngularVelocity
        (
            const setFreePulleyAngularVelocity&
        );

        //- Disallow default bitwise assignment
        void operator=(const setFreePulleyAngularVelocity&);


public:

    //- Runtime type information
    TypeName("setFreePulleyAngularVelocity");


    // Constructors

        //- Construct from components
        setFreePulleyAngularVelocity
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
