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
    myFollowerForceBeamDisplacementNRFvPatchVectorField

Description
	Specifies force perpendicular to the position curve of the beam. Follower load boundary	   		condition

SourceFiles
    followerForceBeamDisplacementNRFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef myFollowerForceBeamDisplacementNRFvPatchVectorField_H
#define myFollowerForceBeamDisplacementNRFvPatchVectorField_H

#include "fvPatchFields.H"
#include "forceBeamDisplacementNRFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class myFollowerForceBeamDisplacementNRFvPatch Declaration
\*---------------------------------------------------------------------------*/

class myFollowerForceBeamDisplacementNRFvPatchVectorField
:
    public forceBeamDisplacementNRFvPatchVectorField
{
    // Private Data

        //- Follower force plane normal
        vectorField followerPlaneNormal_;

        //- Current follower force
        scalarField followerForce_;

        //- Current follower force direction
        vectorField followerForceDir_;
    
        //- Current out of plane force
        scalarField outOfPlaneForce_;

        //- Follower force time series
        interpolationTable<scalar> followerForceSeries_;

        //- Out of plane force time series
        interpolationTable<scalar> outOfPlaneForceSeries_;

        //- Current time index
        label timeIndex_;
    
public:

    //- Runtime type information
    TypeName("myFollowerForceBeamDisplacementNR");

    // Constructors

        //- Construct from patch and internal field
        myFollowerForceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myFollowerForceBeamDisplacementNRFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  myFollowerForceBeamDisplacementNRFvPatchVectorField onto a new patch
        myFollowerForceBeamDisplacementNRFvPatchVectorField
        (
            const myFollowerForceBeamDisplacementNRFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myFollowerForceBeamDisplacementNRFvPatchVectorField
        (
            const myFollowerForceBeamDisplacementNRFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new myFollowerForceBeamDisplacementNRFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myFollowerForceBeamDisplacementNRFvPatchVectorField
        (
            const myFollowerForceBeamDisplacementNRFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new myFollowerForceBeamDisplacementNRFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // //- Evaluate the patch field, sets Updated to false
        // virtual void evaluate
        // (
        //     const Pstream::commsTypes commsType = Pstream::blocking
        // );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
