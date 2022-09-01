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

#include "fvMeshCsrAddressing.H"
#include "labelPair.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshCsrAddressing::fvMeshCsrAddressing(const fvMesh& mesh)
:
    lduAddr_(mesh.lduAddr()),
    cellIndex_
    (
        IOobject
        (
            "cellIndex",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    procNCells_(),
    procCells_(0),
    globalNCellsOffset_(0),
    rowPointers_(),
    columnIndices_(),
    symmRowPointers_(),
    symmColumnIndices_()
{
    label nCells = lduAddr_.size();
    
    if (Pstream::parRun())
    {
        // Set cell index field
        for (label cellI=0; cellI<nCells; cellI++)
        {
            cellIndex_[cellI] = cellI;
        }
        cellIndex_.correctBoundaryConditions();
    
        // Set number of cells per processors
        procNCells_.setSize(Pstream::nProcs());
        procNCells_[Pstream::myProcNo()] = nCells;
        Pstream::gatherList(procNCells_);
        Pstream::scatterList(procNCells_);

        // Set global number of cells offset
        for(label pI=0; pI<Pstream::myProcNo(); pI++)
        {
            globalNCellsOffset_ += procNCells_[pI];
        }

        // Find cells next to the processor boundary
        // Note: assumed that only one processor face per cell exists
        procCells_.setSize(nCells, FixedList<label, 3>(-1));
        forAll(mesh.boundary(), patchI)
        {
            if (isA<processorFvPatch>(mesh.boundary()[patchI]))
            {
                const processorFvPatch& curProcPatch =
                    refCast<const processorFvPatch>(mesh.boundary()[patchI]);

                const unallocLabelList& fc =
                    mesh.boundary()[patchI].faceCells();

                label neiProcNo = curProcPatch.neighbProcNo();

                // Determina global number of cells for previous procesors
                label curGlobalNCellsOffset = 0;
                for(label pI=0; pI<neiProcNo; pI++)
                {
                    curGlobalNCellsOffset += procNCells_[pI];
                }

                scalarField neiLocalCellIndices =
                    cellIndex_.boundaryField()[patchI].patchNeighbourField();

                scalarField neiGlobalCellIndices =
                    curGlobalNCellsOffset + neiLocalCellIndices;

                forAll(fc, faceI)
                {
                    procCells_[fc[faceI]][0] = patchI;
                    procCells_[fc[faceI]][1] = faceI;
                    procCells_[fc[faceI]][2] = neiGlobalCellIndices[faceI];
                }
            }
        }
    }
    else
    {
        procNCells_.setSize(1);
        procNCells_[0] = nCells;
    }

    // Unsymmetric matrix addressing
    rowPointers_.setSize(nCells+1, -1);
    DynamicList<label> columns;

    rowPointers_[0] = 0;
    for (label rowI=1; rowI<rowPointers_.size(); rowI++)
    {
        // Add lower nei processor cells
        if (procCells_.size())
        {
            if
            (
                procCells_[rowI-1][0] != -1
             && procCells_[rowI-1][2] < globalNCellsOffset_
            )
            {
                columns.append(procCells_[rowI-1][2]);
            }
        }

        // Lower
        for
        (
            label j=lduAddr_.losortStartAddr()[rowI-1];
            j<lduAddr_.losortStartAddr()[rowI];
            j++
        )
        {
            columns.append
            (
                lduAddr_.lowerAddr()[lduAddr_.losortAddr()[j]]
              + globalNCellsOffset_
            );
        }

        // Diagonal
        columns.append((rowI-1) + globalNCellsOffset_);

        // Upper
        for
        (
            label j=lduAddr_.ownerStartAddr()[rowI-1];
            j<lduAddr_.ownerStartAddr()[rowI];
            j++
        )
        {
            columns.append(lduAddr_.upperAddr()[j] + globalNCellsOffset_);
        }

        // Add upper nei processor cells
        if (procCells_.size())
        {
            if
            (
                procCells_[rowI-1][0] != -1
             && procCells_[rowI-1][2] >= (globalNCellsOffset_ + nCells)
            )
            {
                columns.append(procCells_[rowI-1][2]);
            }
        }

        rowPointers_[rowI] = columns.size();
    }
    columnIndices_ = columns;

    // Info << lduAddr_.ownerStartAddr() << endl;
    // Info << lduAddr_.upperAddr() << endl;
    
    // Info << lduAddr_.losortStartAddr() << endl;
    // Info << lduAddr_.losortAddr() << endl;

    // Info << rowPointers_ << endl;
    // Info << columnIndices_ << endl;

    // // Symmetric matrix addressing
    // symmRowPointers_.setSize(lduAddr_.size()+1, -1);
    // columns.clear();

    // symmRowPointers_[0] = 0;
    // for (label rowI=1; rowI<rowPointers_.size(); rowI++)
    // {
    //     // Diagonal
    //     columns.append(rowI-1);

    //     // Upper
    //     for
    //     (
    //         label j=lduAddr_.ownerStartAddr()[rowI-1];
    //         j<lduAddr_.ownerStartAddr()[rowI];
    //         j++
    //     )
    //     {
    //         columns.append(lduAddr_.upperAddr()[j]);
    //     }

    //     symmRowPointers_[rowI] = columns.size();
    // }
    // symmColumnIndices_ = columns;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
