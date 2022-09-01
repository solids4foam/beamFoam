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
    Foam::fvMeshCsrAddressing

Description
    Foam::fvMeshCsrAddressing

SourceFiles
    fvMeshCsrAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshCsrAddressing_H
#define fvMeshCsrAddressing_H

#include "lduAddressing.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fvMeshCsrAddressing Declaration
\*---------------------------------------------------------------------------*/

class fvMeshCsrAddressing
{
    // Private data

        //- Reference to ldu addressing
        const lduAddressing& lduAddr_;

        //- Cell index field
        volScalarField cellIndex_;

        //- Number of cells per processors
        labelList procNCells_;

        //- Cells next to processor boundary
        List<FixedList<label, 3> > procCells_;

        //- Global number of cells offset
        label globalNCellsOffset_;
    
        //- Row pointers
        labelList rowPointers_;

        //- Column indices
        labelList columnIndices_;

        //- Row pointers for symmetric matrix
        labelList symmRowPointers_;

        //- Column indices for symmetric matrix
        labelList symmColumnIndices_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvMeshCsrAddressing(const fvMeshCsrAddressing&);

        //- Disallow default bitwise assignment
        void operator=(const fvMeshCsrAddressing&);

public:

    // Constructors

        //- Construct from components
        fvMeshCsrAddressing(const fvMesh& mesh);

        //- Destructor
        virtual ~fvMeshCsrAddressing()
        {}

    // Member Functions

        //- Return cell index field
        const volScalarField& cellIndex() const
        {
            return cellIndex_;
        }

        //- Number of cells per processors
        const labelList& procNCells() const
        {
            return procNCells_;
        }

        //- Cells next to processor boundary
        const List<FixedList<label, 3> >& procCells() const
        {
            return procCells_;
        }

        //- Global number of cells offset
        const label& globalNCellsOffset() const
        {
            return globalNCellsOffset_;
        }
    
        //- Return row pointers
        const labelList& rowPointers() const
        {
            return rowPointers_;
        }

        //- Return column indices
        const labelList& columnIndices() const
        {
            return columnIndices_;
        }
    
        //- Return row pointers
        const labelList& symmRwPointers() const
        {
            return symmRowPointers_;
        }

        //- Return column indices
        const labelList& symmColumnIndices() const
        {
            return symmColumnIndices_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
