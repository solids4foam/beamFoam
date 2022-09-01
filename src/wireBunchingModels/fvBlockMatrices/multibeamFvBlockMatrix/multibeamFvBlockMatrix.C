/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "multibeamFvBlockMatrix.H"
#include "surfaceInterpolationScheme.H"
#include "beamModel.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multibeamFvBlockMatrix, 0);
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

Foam::label Foam::multibeamFvBlockMatrix::
lowerPointContact(const label bI, const label segI) const
{
    label pointContactIndex = -1;

    forAll(beam_.contact().pointContacts(), pcI)
    {
        if
        (
            beam_.contact().pointContacts()[pcI].firstBeam() == bI
         && beam_.contact().pointContacts()[pcI].firstBeamSegment() == segI
        )
        {
            pointContactIndex = pcI;
            break;
        }
    }

    return pointContactIndex;
}

Foam::label Foam::multibeamFvBlockMatrix::
upperPointContact(const label bI, const label segI) const
{
    label pointContactIndex = -1;
    
    forAll(beam_.contact().pointContacts(), pcI)
    {
        if
        (
            beam_.contact().pointContacts()[pcI].secondBeam() == bI
         && beam_.contact().pointContacts()[pcI].secondBeamSegment() == segI
        )
        {
            pointContactIndex = pcI;
            break;
        }
    }

    return pointContactIndex;
}


void Foam::multibeamFvBlockMatrix::makeCsrMatrix() const
{
    if (debug)
    {
        Info<< "multibeamFvBlockMatrix::makeCsrMatrix() const : "
            << "Covertig block CSR to CSR matrix"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        rowPointersPtr_
     || columnIndicesPtr_
     || coeffsPtr_
    )
    {
        FatalErrorIn
        (
            "multibeamFvBlockMatrix::makeCsrMatrix() const"
        )
            << "CSR matrix already exists"
            << abort(FatalError);
    }

    const lduAddressing& lduAddr = this->lduAddr();

    // Block size
    label blockSize = 6;

    // Number of cells (eqs)
    label nCells = lduAddr.size();

    // Block CSR addressing
    const labelList& brp = blockRowPointers();
    const labelList& bci = blockColumnIndices();
    const scalarField& bCoeffs = blockCoeffs();

    // Lagrange multipliers part
    label nLM = 0;
    labelList LMtoc = activeConicalPulleyContacts().toc();
    if (LMtoc.size())
    {
        // Double contact point for conical pulley
        nLM += 2*LMtoc.size(); 
    }

    // CSR row pointers
    label nRows = nCells*blockSize + nLM;
    rowPointersPtr_ = new labelList(nRows+1, -1);
    labelList& rp = *rowPointersPtr_;

    // CSR coefficients
    label nCoeffs = bCoeffs.size() + nLM*(3+4); 
    coeffsPtr_ = new scalarField(nCoeffs, 0);
    scalarField& coeffs = *coeffsPtr_;

    // CSR column indices
    columnIndicesPtr_ = new labelList(coeffs.size(), 0);
    labelList& ci = *columnIndicesPtr_;

    rp[0] = 0;
    label k = 0;
    label i = 0;
    for (label rowI=1; rowI<brp.size(); rowI++)
    {
        label colStart = brp[rowI-1];
        label colEnd = brp[rowI];

        label glSegI = rowI-1;
        label LMindex = findIndex(LMtoc, glSegI);
        vectorField contactForceDirection(2, vector::zero);

        if (nLM)
        {
            if (LMindex != -1)
            {
                const FixedList<label,3>& curContact =
                    activeConicalPulleyContacts()[glSegI];
                label bI = curContact[0];
                label segI = curContact[1];
                label pulleyI = curContact[2];

                contactForceDirection =
                    beam_.contact().conicalPulleyContacts()[bI][segI][pulleyI]
                   .normalDirection();
            }
        }

        for (label rI=0; rI<blockSize; rI++)
        {
            for (label colI=colStart; colI<colEnd; colI++)
            {
                label j = blockSize*bci[colI];
                label c = colI*(blockSize*blockSize) + rI*blockSize;

                for (label cI=0; cI<blockSize; cI++)
                {
                    ci[k] = j++;
                    coeffs[k++] = bCoeffs[c++];
                }
            }

            if (nLM)
            {
                if (LMindex != -1)
                {
                    if (rI == 0)
                    {
                        for (label cpI=0; cpI<2; cpI++)
                        {
                            ci[k] = nCells*blockSize + (2*LMindex + cpI);
                            coeffs[k++] =
                                contactForceDirection[cpI].x()*
                                beam_.L()[glSegI];
                        }
                    }
                    else if (rI == 1)
                    {
                        for (label cpI=0; cpI<2; cpI++)
                        {
                            ci[k] = nCells*blockSize + (2*LMindex + cpI);
                            coeffs[k++] =
                                contactForceDirection[cpI].y()*
                                beam_.L()[glSegI];
                        }
                    }
                    else if (rI == 2)
                    {
                        for (label cpI=0; cpI<2; cpI++)
                        {
                            ci[k] = nCells*blockSize + (2*LMindex + cpI);
                            coeffs[k++] =
                                contactForceDirection[cpI].z()*
                                beam_.L()[glSegI];
                        }
                    }
                }
            }

            rp[++i] = k;
        }
    }

    // Add lagrange multipliers
    if (nLM)
    {
        label LMindex = 0;
        (*rhsPtr_).setSize(rhs().size() + LMtoc.size()*2);
        (*solutionPtr_).setSize(rhs().size() + LMtoc.size()*2);
        for
        (
            Map<FixedList<label,3> >::const_iterator iter =
                activeConicalPulleyContacts().begin();
            iter != activeConicalPulleyContacts().end();
            ++iter
        )
        {
            label glSegI = LMtoc[LMindex];

            label bI = (*iter)[0];
            label segI = (*iter)[1];
            label pulleyI = (*iter)[2];

            const vectorField& n =
                beam_.contact().conicalPulleyContacts()[bI][segI][pulleyI]
               .normalDirection();

            scalarField fn =
                mag
                (
                    beam_.contact().conicalPulleyContacts()[bI][segI][pulleyI]
                   .normalContactForce()
                );

            const scalarField& gn =
                beam_.contact().conicalPulleyContacts()[bI][segI][pulleyI]
               .normalGap();

            // Info << bI << ", " << segI << ", " << pulleyI << endl;
            // Info << "gn: " << gn << endl; 
            // Info << "fn: " << fn << endl;
            // Info << "n: " << n << endl;
       
            for (label cp=0; cp<2; cp++)
            {
                ci[k] = glSegI*blockSize+0;
                coeffs[k++] = n[cp].x();
                // coeffs[k++] = fn[cp]*n[cp].x();
                
                ci[k] = glSegI*blockSize+1;
                coeffs[k++] = n[cp].y();
                // coeffs[k++] = fn[cp]*n[cp].y();
            
                ci[k] = glSegI*blockSize+2;
                coeffs[k++] = n[cp].z();
                // coeffs[k++] = fn[cp]*n[cp].z();
                
                // ci[k] = nCells*blockSize + (2*LMindex + cp);
                // coeffs[k++] = sqr(gn[cp])/fn[cp];

                (*rhsPtr_)[i] = -gn[cp];
                // (*rhsPtr_)[i] = -gn[cp]*fn[cp];
                (*solutionPtr_)[i] = 0;
                
                rp[++i] = k;
            }

            LMindex++;
        }
    }
}


void Foam::multibeamFvBlockMatrix::makeBlockCsrMatrix() const
{
    if (debug)
    {
        Info<< "multibeamFvBlockMatrix::makeBlockCsrMatrix() const : "
            << "Covertig LDU to CSR matrix"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        blockRowPointersPtr_
     || blockColumnIndicesPtr_
     || blockCoeffsPtr_
     || rhsPtr_
     || solutionPtr_
    )
    {
        FatalErrorIn
        (
            "multibeamFvBlockMatrix::makeCsrMatrix() const"
        )
            << "Block CSR matrix already exists"
            << abort(FatalError);
    }

    // Matrix addressing
    const lduAddressing& lduAddr = this->lduAddr();
    const unallocLabelList& losortAddr = lduAddr.losortAddr();
    const unallocLabelList& losortStartAddr = lduAddr.losortStartAddr();
    const unallocLabelList& ownerStartAddr = lduAddr.ownerStartAddr();

    // Grab matrix diagonal and off-diagonals
    const Field<tensor6>& d = this->diag().asSquare();
    const Field<tensor6>& l = this->lower().asSquare();
    const Field<tensor6>& u = this->upper().asSquare();
    const Field<vector6>& blockRhs = this->source();

    // Block size
    label blockSize = 6;

    // Number of cells (eqs)
    label nCells = lduAddr.size();

    // Number of beams
    label nBeams = beam_.mesh().cellZones().size();

    // CSR addressing without inter-beam coupling
    const labelList& oldBlockRowPointers = beam_.csrAddr().rowPointers();
    const labelList& oldBlockColumnIndices = beam_.csrAddr().columnIndices();

    // CSR row pointers
    blockRowPointersPtr_ = new labelList(nCells+1, -1);
    labelList& rp = *blockRowPointersPtr_;

    // CSR column indices
    DynamicList<label> columns;

    rp[0] = 0;
    for (label rowI=1; rowI<rp.size(); rowI++)
    {
        // Add point contact coupling addressing
        // Assuming that line and point contacts can not occur together
        if (pointContactCoeffs_.size() && !Pstream::parRun())
        {            
            label bI = beam_.whichBeam(rowI-1);
            label segI = beam_.whichSegment(rowI-1);

            label pcI = upperPointContact(bI, segI);

            if (pcI != -1)
            {
                const label lowerBeam = 
                    beam_.contact().pointContacts()[pcI]
                   .firstBeam();
                
                const label lowerSeg =
                    beam_.contact().pointContacts()[pcI]
                   .firstBeamLowerSegment();
                
                const label upperSeg =
                    beam_.contact().pointContacts()[pcI]
                   .firstBeamUpperSegment();

                if (lowerSeg != upperSeg)
                {
                    const label lowerCell =
                        beam_.whichCell(lowerBeam, lowerSeg);
                    const label upperCell =
                        beam_.whichCell(lowerBeam, upperSeg);

                    columns.append(lowerCell);
                    columns.append(upperCell);
                }
                else
                {
                    const label lowerCell =
                        beam_.whichCell(lowerBeam, lowerSeg);

                    columns.append(lowerCell);
                }
            }
        }

        // Add lower beams coupling addressing (for line contact only)
        if (beam_.beamContactActive())
        {
            label localCellIndex = rowI-1;
            label gcI = beam_.globalCellIndex(localCellIndex);
            
            // label globalCellIndex = localCellIndex
            //   + beam_.csrAddr().globalNCellsOffset();

            label bI = beam_.whichBeam(gcI);
            label segI = beam_.whichSegment(gcI);
            // label bI = beam_.whichBeam(rowI-1);
            // label segI = beam_.whichSegment(rowI-1);

            const lineContactList& curLineContacts =
	        beam_.contact().lineContacts()[bI][segI];

            for (label nbI=0; nbI<bI; nbI++)
            {
                if (mag(curLineContacts[nbI].normalContactForce()) > SMALL)
                {
                    const label neiLowerSegI =
                        curLineContacts[nbI].secondBeamLowerSegment();
                    const label neiUpperSegI =
                        curLineContacts[nbI].secondBeamUpperSegment();

                    const label neiLowerCell =
                        beam_.whichCell(nbI, neiLowerSegI);
                    const label locNeiLowerCell =
                        beam_.procLocalCellIndex(neiLowerCell).second();

                    const label neiUpperCell =
                        beam_.whichCell(nbI, neiUpperSegI);
                    const label locNeiUpperCell =
                        beam_.procLocalCellIndex(neiUpperCell).second();

                    // Pout << "lower " << localCellIndex << ", "
                    //      << locNeiLowerCell << ", "
                    //      << locNeiUpperCell << endl;
                        
                    if (neiLowerSegI != neiUpperSegI)
                    {
                        columns.append(locNeiLowerCell);
                        columns.append(locNeiUpperCell);
                    }
                    else
                    {
                        columns.append(locNeiLowerCell);
                    }
                    
                    // if (neiLowerSegI != neiUpperSegI)
                    // {
                    //     const label neiLowerCell =
                    //         beam_.whichCell(nbI, neiLowerSegI);
                    //     const label neiUpperCell =
                    //         beam_.whichCell(nbI, neiUpperSegI);

                    //     columns.append(neiLowerCell);
                    //     columns.append(neiUpperCell);
                    // }
                    // else
                    // {
                    //     const label neiLowerCell =
                    //         beam_.whichCell(nbI, neiLowerSegI);

                    //     columns.append(neiLowerCell);
                    // }
                }
            }
        }

        // Add old columns
        for
        (
            label cI = oldBlockRowPointers[rowI-1];
            cI<oldBlockRowPointers[rowI];
            cI++
        )
        {
            columns.append(oldBlockColumnIndices[cI]);
        }

        // Add upper beams coupling (line contact only)
        if (beam_.beamContactActive())
        {
            label localCellIndex = rowI-1;
            label gcI = beam_.globalCellIndex(localCellIndex);

            // label globalCellIndex = localCellIndex
            //   + beam_.csrAddr().globalNCellsOffset();

            label bI = beam_.whichBeam(gcI);
            label segI = beam_.whichSegment(gcI);
            // label bI = beam_.whichBeam(rowI-1);
            // label segI = beam_.whichSegment(rowI-1);

            const lineContactList& curLineContacts =
	        beam_.contact().lineContacts()[bI][segI];

            for (label nbI=bI+1; nbI<nBeams; nbI++)
            {
                if (mag(curLineContacts[nbI].normalContactForce()) > SMALL)
                {                    
                    const label neiLowerSegI =
                        curLineContacts[nbI].secondBeamLowerSegment();
                    const label neiUpperSegI =
                        curLineContacts[nbI].secondBeamUpperSegment();

                    const label neiLowerCell =
                        beam_.whichCell(nbI, neiLowerSegI);
                    const label locNeiLowerCell =
                        beam_.procLocalCellIndex(neiLowerCell).second();

                    const label neiUpperCell =
                        beam_.whichCell(nbI, neiUpperSegI);
                    const label locNeiUpperCell =
                        beam_.procLocalCellIndex(neiUpperCell).second();
                        
                    // Pout << "upper " << localCellIndex << ", "
                    //      << locNeiLowerCell << ", "
                    //      << locNeiUpperCell << endl;
                    
                    if (neiLowerSegI != neiUpperSegI)
                    {
                        columns.append(locNeiLowerCell);
                        columns.append(locNeiUpperCell);
                    }
                    else
                    {
                        columns.append(locNeiLowerCell);
                    }
                }
            }
        }

        // Add point contact coupling addressing
        // Assuming that line and point contacts can not occur together
        if (pointContactCoeffs_.size() && !Pstream::parRun())
        {            
            label bI = beam_.whichBeam(rowI-1);
            label segI = beam_.whichSegment(rowI-1);

            label pcI = lowerPointContact(bI, segI);

            if (pcI != -1)
            {
                const label upperBeam = 
                    beam_.contact().pointContacts()[pcI]
                   .secondBeam();

                const label lowerSeg =
                    beam_.contact().pointContacts()[pcI]
                   .secondBeamLowerSegment();

                const label upperSeg =
                    beam_.contact().pointContacts()[pcI]
                   .secondBeamUpperSegment();

                if (lowerSeg != upperSeg)
                {
                    const label lowerCell =
                        beam_.whichCell(upperBeam, lowerSeg);
                    const label upperCell =
                        beam_.whichCell(upperBeam, upperSeg);

                    columns.append(lowerCell);
                    columns.append(upperCell);
                }
                else
                {
                    const label lowerCell =
                        beam_.whichCell(upperBeam, lowerSeg);

                    columns.append(lowerCell);
                }
            }
        }

        rp[rowI] = columns.size();

        // if (Pstream::myProcNo() == 0)
        // {
        //     if (rowI-1 == 4)
        //     {
        //         const label jStart = rp[rowI-1];
        //         const label jEnd = rp[rowI];

        //         Pout << rowI-1 << endl;
        //         for (label j=jStart; j<jEnd; j++)
        //         {
        //             const label cI = columns[j];

        //             Pout << cI << " ";
        //         }
        //         Pout << endl;
        //     }
        // }
    }

    blockColumnIndicesPtr_ = new labelList(columns);
    const labelList& ci = *blockColumnIndicesPtr_;

    // Calculate coeffs
    blockCoeffsPtr_ = new scalarField(blockSize*blockSize*ci.size(), 0);
    scalarField& blockCoeffs = *blockCoeffsPtr_;
    label coeffI = 0;
    for (label rowI=1; rowI<rp.size(); rowI++)
    {
        bool activeContact = false;
        
        // Add point contact coupling coefficients
        // Assuming that line and point contacts can not occur together
        if (pointContactCoeffs_.size() && !Pstream::parRun())
        {            
            label bI = beam_.whichBeam(rowI-1);
            label segI = beam_.whichSegment(rowI-1);

            label pcI = upperPointContact(bI, segI);

            if (pcI != -1)
            {
                const label lowerSeg =
                    beam_.contact().pointContacts()[pcI]
                   .firstBeamLowerSegment();

                const label upperSeg =
                    beam_.contact().pointContacts()[pcI]
                   .firstBeamUpperSegment();

                if (lowerSeg != upperSeg)
                {
                    const tensor6& curLowerCoeff =
                        pointContactCoeffs_[pcI][1].first();

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }
                    
                    const tensor6& curUpperCoeff =
                        pointContactCoeffs_[pcI][1].second();

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curUpperCoeff(ii, jj);
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curUpperCoeff(ii, jj);
                            }
                        }
                    }
                }
                else
                {
                    const tensor6& curLowerCoeff =
                        pointContactCoeffs_[pcI][1].first();

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }
                }
            }
        }

        // Add lower beams coupling coeffs (line contact)
        if (beam_.beamContactActive())
        {
            label localCellIndex = rowI-1;
            label gcI = beam_.globalCellIndex(localCellIndex);

            // label globalCellIndex = localCellIndex
            //   + beam_.csrAddr().globalNCellsOffset();

            label bI = beam_.whichBeam(gcI);
            label segI = beam_.whichSegment(gcI);

            // label bI = beam_.whichBeam(rowI-1);
            // label segI = beam_.whichSegment(rowI-1);

            const lineContactList& curLineContacts =
	        beam_.contact().lineContacts()[bI][segI];

            for (label nbI=0; nbI<bI; nbI++)
            {
                if (mag(curLineContacts[nbI].normalContactForce()) > SMALL)
                {
                    activeContact = true;
                    
                    const label neiLowerSegI =
                        curLineContacts[nbI].secondBeamLowerSegment();
                    const label neiUpperSegI =
                        curLineContacts[nbI].secondBeamUpperSegment();

                    if (neiLowerSegI != neiUpperSegI)
                    {
                        const tensor6& curLowerCoeff =
                            lineContactCoeffs_[bI][segI][nbI].first();

                        if (rowMajor_)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                for (label jj=0; jj<6; jj++)
                                {
                                    // row-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }
                        else
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                for (label ii=0; ii<6; ii++)
                                {
                                    // column-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }

                        const tensor6& curUpperCoeff =
                            lineContactCoeffs_[bI][segI][nbI].second();

                        if (rowMajor_)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                for (label jj=0; jj<6; jj++)
                                {
                                    // row-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curUpperCoeff(ii,jj);
                                }
                            }
                        }
                        else
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                for (label ii=0; ii<6; ii++)
                                {
                                    // column-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curUpperCoeff(ii,jj);
                                }
                            }
                        }
                    }
                    else
                    {
                        const tensor6& curLowerCoeff =
                            lineContactCoeffs_[bI][segI][nbI].first();

                        if (rowMajor_)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                for (label jj=0; jj<6; jj++)
                                {
                                    // row-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }
                        else
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                for (label ii=0; ii<6; ii++)
                                {
                                    // column-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Add lower nei processor cells
        if (beam_.csrAddr().procCells().size())
        {
            if
            (
                beam_.csrAddr().procCells()[rowI-1][0] != -1
             && beam_.csrAddr().procCells()[rowI-1][2]
              < beam_.csrAddr().globalNCellsOffset()
            )
            {
                label patchI = beam_.csrAddr().procCells()[rowI-1][0];
                label faceI = beam_.csrAddr().procCells()[rowI-1][1];
                
                const tensor6& curCoeff = 
                    this->coupleUpper()[patchI].asSquare()[faceI];
                
                if (rowMajor_)
                {
                    for (label ii=0; ii<6; ii++)
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            // row-major order for 0-based indexing
                            blockCoeffs[coeffI++] = -curCoeff(ii,jj);
                        }
                    }
                }
                else
                {
                    for (label jj=0; jj<6; jj++)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            // column-major order for 0-based indexing
                            blockCoeffs[coeffI++] = -curCoeff(ii,jj);
                        }
                    }
                }
            }
        }

        // Lower
        for
        (
            label j=losortStartAddr[rowI-1];
            j<losortStartAddr[rowI];
            j++
        )
        {
            label faceI = losortAddr[j];
            
            const tensor6& curCoeff = l[faceI];

            if (rowMajor_)
            {
                for (label ii=0; ii<6; ii++)
                {
                    for (label jj=0; jj<6; jj++)
                    {
                        // row-major order for 0-based indexing
                        blockCoeffs[coeffI++] = curCoeff(ii,jj);
                    }
                }
            }
            else
            {
                for (label jj=0; jj<6; jj++)
                {
                    for (label ii=0; ii<6; ii++)
                    {
                        // column-major order for 0-based indexing
                        blockCoeffs[coeffI++] = curCoeff(ii,jj);
                    }
                }
            }
        }        

        // Diagonal
        {
            const tensor6& curCoeff = d[rowI-1];

            if (rowMajor_)
            {
                for (label ii=0; ii<6; ii++)
                {
                    for (label jj=0; jj<6; jj++)
                    {
                        // row-major order for 0-based indexing
                        blockCoeffs[coeffI++] = curCoeff(ii,jj);
                    }
                }
            }
            else
            {
                for (label jj=0; jj<6; jj++)
                {
                    for (label ii=0; ii<6; ii++)
                    {
                        // column-major order for 0-based indexing
                        blockCoeffs[coeffI++] = curCoeff(ii,jj);
                    }
                }
            }
        }

        // Upper
        for
        (
            label j=ownerStartAddr[rowI-1];
            j<ownerStartAddr[rowI];
            j++
        )
        {
            label faceI = j;
            
            const tensor6& curCoeff = u[faceI];

            if (rowMajor_)
            {
                for (label ii=0; ii<6; ii++)
                {
                    for (label jj=0; jj<6; jj++)
                    {
                        // row-major order for 0-based indexing
                        blockCoeffs[coeffI++] = curCoeff(ii,jj);
                    }
                }
            }
            else
            {
                for (label jj=0; jj<6; jj++)
                {
                    for (label ii=0; ii<6; ii++)
                    {
                        // column-major order for 0-based indexing
                        blockCoeffs[coeffI++] = curCoeff(ii,jj);
                    }
                }
            }
        }

        // Add upper nei processor cells
        if (beam_.csrAddr().procCells().size())
        {
            if
            (
                beam_.csrAddr().procCells()[rowI-1][0] != -1
             && beam_.csrAddr().procCells()[rowI-1][2]
             >= (beam_.csrAddr().globalNCellsOffset() + nCells)
            )
            {
                label patchI = beam_.csrAddr().procCells()[rowI-1][0];
                label faceI = beam_.csrAddr().procCells()[rowI-1][1];

                const tensor6& curCoeff = 
                    this->coupleUpper()[patchI].asSquare()[faceI];

                if (rowMajor_)
                {
                    for (label ii=0; ii<6; ii++)
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            // row-major order for 0-based indexing
                            blockCoeffs[coeffI++] = -curCoeff(ii,jj);
                        }
                    }
                }
                else
                {
                    for (label jj=0; jj<6; jj++)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            // column-major order for 0-based indexing
                            blockCoeffs[coeffI++] = -curCoeff(ii,jj);
                        }
                    }
                }
            }
        }

        // Add upper beams coupling coeffs
        if (beam_.beamContactActive())
        {
            label localCellIndex = rowI-1;
            label gcI = beam_.globalCellIndex(localCellIndex);

            // label globalCellIndex = localCellIndex
            //   + beam_.csrAddr().globalNCellsOffset();

            label bI = beam_.whichBeam(gcI);
            label segI = beam_.whichSegment(gcI);

            // label bI = beam_.whichBeam(rowI-1);
            // label segI = beam_.whichSegment(rowI-1);

            const lineContactList& curLineContacts =
	        beam_.contact().lineContacts()[bI][segI];

            for (label nbI=bI+1; nbI<nBeams; nbI++)
            {
                if (mag(curLineContacts[nbI].normalContactForce()) > SMALL)
                {
                    activeContact = true;
                    
                    const label neiLowerSegI =
                        curLineContacts[nbI].secondBeamLowerSegment();
                    const label neiUpperSegI =
                        curLineContacts[nbI].secondBeamUpperSegment();

                    if (neiLowerSegI != neiUpperSegI)
                    {
                        const tensor6& curLowerCoeff =
                            lineContactCoeffs_[bI][segI][nbI].first();

                        if (rowMajor_)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                for (label jj=0; jj<6; jj++)
                                {
                                    // row-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }
                        else
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                for (label ii=0; ii<6; ii++)
                                {
                                    // column-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }

                        const tensor6& curUpperCoeff =
                            lineContactCoeffs_[bI][segI][nbI].second();

                        if (rowMajor_)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                for (label jj=0; jj<6; jj++)
                                {
                                    // row-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curUpperCoeff(ii,jj);
                                }
                            }
                        }
                        else
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                for (label ii=0; ii<6; ii++)
                                {
                                    // column-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curUpperCoeff(ii,jj);
                                }
                            }
                        }
                    }
                    else
                    {
                        const tensor6& curLowerCoeff =
                            lineContactCoeffs_[bI][segI][nbI].first();

                        if (rowMajor_)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                for (label jj=0; jj<6; jj++)
                                {
                                    // row-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }
                        else
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                for (label ii=0; ii<6; ii++)
                                {
                                    // column-major order for 0-based indexing
                                    blockCoeffs[coeffI++] = curLowerCoeff(ii,jj);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Add point contact coupling coefficients
        // Assuming that line and point contacts can not occur together
        if (pointContactCoeffs_.size() && !Pstream::parRun())
        {            
            label bI = beam_.whichBeam(rowI-1);
            label segI = beam_.whichSegment(rowI-1);

            label pcI = lowerPointContact(bI, segI);

            if (pcI != -1)
            {
                const label lowerSeg =
                    beam_.contact().pointContacts()[pcI]
                   .secondBeamLowerSegment();

                const label upperSeg =
                    beam_.contact().pointContacts()[pcI]
                   .secondBeamUpperSegment();

                if (lowerSeg != upperSeg)
                {
                    const tensor6& curLowerCoeff =
                        pointContactCoeffs_[pcI][0].first();

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }

                    const tensor6& curUpperCoeff =
                        pointContactCoeffs_[pcI][0].second();

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curUpperCoeff(ii, jj);
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curUpperCoeff(ii, jj);
                            }
                        }
                    }
                }
                else
                {
                    const tensor6& curLowerCoeff =
                        pointContactCoeffs_[pcI][0].first();

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                blockCoeffs[coeffI++] = curLowerCoeff(ii, jj);
                            }
                        }
                    }
                }
            }
        }

        // Sorting if necessary
        if (beam_.csrAddr().procCells().size())
        {
            if
            (
                (beam_.csrAddr().procCells()[rowI-1][0] != -1)  // proc patch
             && activeContact
            )
            {
                const label jStart = rp[rowI-1];
                const label jEnd = rp[rowI];

                label nTerms = jEnd - jStart;

                SortableList<label> curCols(nTerms);
                Field<tensor6> curCoeffs(nTerms, tensor6::zero);

                label bs = 6*6;

                label colI=0;
                for (label j=jStart; j<jEnd; j++)
                {
                    curCols[colI] = ci[j];

                    label jBlock = j*bs;

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                curCoeffs[colI](ii, jj) = blockCoeffs[jBlock++];
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                curCoeffs[colI](ii, jj) = blockCoeffs[jBlock++];
                            }
                        }
                    }
                    
                    colI++;
                }

                curCols.sort();

                colI=0;
                for (label j=jStart; j < jEnd; j++)
                {
                    (*blockColumnIndicesPtr_)[j] = curCols[colI];

                    label jBlock = j*bs;

                    if (rowMajor_)
                    {
                        for (label ii=0; ii<6; ii++)
                        {
                            for (label jj=0; jj<6; jj++)
                            {
                                // row-major order for 0-based indexing
                                blockCoeffs[jBlock++] =
                                    curCoeffs[curCols.indices()[colI]](ii, jj);
                            }
                        }
                    }
                    else
                    {
                        for (label jj=0; jj<6; jj++)
                        {
                            for (label ii=0; ii<6; ii++)
                            {
                                // column-major order for 0-based indexing
                                blockCoeffs[jBlock++] =
                                    curCoeffs[curCols.indices()[colI]](ii, jj);
                            }
                        }
                    }

                    colI++;
                }

                // Pout << (rowI-1) << endl;
                // Pout << curCols << endl;
            }
        }
    }

    // Calculate rhs and solutin vector
    rhsPtr_ = new scalarField(blockSize*blockRhs.size());
    scalarField& rhs = *rhsPtr_;

    solutionPtr_ = new scalarField(blockSize*blockRhs.size());
    scalarField& solution = *solutionPtr_;

    const Field<vector6>& blockSolution = this->psi().internalField();
    label k = 0;
    for (label i=0; i<nCells; i++)
    {
        for (label j=0; j<blockSize; j++)
        {
            rhs[k] = blockRhs[i](j);
            solution[k] = blockSolution[i](j);
            k++;
        }
    }

    // Pout << "done" << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::multibeamFvBlockMatrix::multibeamFvBlockMatrix
(
    GeometricField<vector6, fvPatchField, volMesh>& psi,
    beamModel& beam
)
:
    fvBlockVector6Matrix(psi),
    beam_(beam),
    lineContactCoeffs_(),
    pointContactCoeffs_(),
    activeConicalPulleyContacts_(),
    blockRowPointersPtr_(NULL),
    blockColumnIndicesPtr_(NULL),
    blockCoeffsPtr_(NULL),
    rowPointersPtr_(NULL),
    columnIndicesPtr_(NULL),
    coeffsPtr_(NULL),
    rhsPtr_(NULL),
    solutionPtr_(NULL),
    rowMajor_(true)
{
    label nBeams = beam_.mesh().cellZones().size();

    if (nBeams > 1)
    {
        // Initialiye line contact off-dagonal coeffs
        lineContactCoeffs_.setSize(nBeams);
        forAll(lineContactCoeffs_, bI)
        {
            lineContactCoeffs_.set
            (
                bI,
                new tensor6PairListList
                (
                    beam_.contact().splines()[bI].nSegments()
                )
            );

            forAll(lineContactCoeffs_[bI], segI)
            {
                lineContactCoeffs_[bI][segI] =
                    tensor6PairList
                    (
                        nBeams,
                        tensor6Pair(tensor6::zero, tensor6::zero)
                    );
            }
        }
    }
}


Foam::multibeamFvBlockMatrix::~multibeamFvBlockMatrix()
{
    // Delete demand driven data
    deleteDemandDrivenData(blockRowPointersPtr_);
    deleteDemandDrivenData(blockColumnIndicesPtr_);
    deleteDemandDrivenData(blockCoeffsPtr_);
    deleteDemandDrivenData(rowPointersPtr_);
    deleteDemandDrivenData(columnIndicesPtr_);
    deleteDemandDrivenData(coeffsPtr_);
    deleteDemandDrivenData(rhsPtr_);
    deleteDemandDrivenData(solutionPtr_);
}


Foam::PtrList<Foam::tensor6PairListList>&
Foam::multibeamFvBlockMatrix::lineContactCoeffs()
{
    // Delete old CSR matrix
    deleteDemandDrivenData(blockRowPointersPtr_);
    deleteDemandDrivenData(blockColumnIndicesPtr_);
    deleteDemandDrivenData(blockCoeffsPtr_);
    deleteDemandDrivenData(rowPointersPtr_);
    deleteDemandDrivenData(columnIndicesPtr_);
    deleteDemandDrivenData(coeffsPtr_);
    deleteDemandDrivenData(rhsPtr_);
    deleteDemandDrivenData(solutionPtr_);

    // Set coeffs to zero
    label nBeams = beam_.mesh().cellZones().size();
    forAll(lineContactCoeffs_, bI)
    {
        forAll(lineContactCoeffs_[bI], segI)
        {
            lineContactCoeffs_[bI][segI] =
                tensor6PairList
                (
                    nBeams,
                    tensor6Pair(tensor6::zero, tensor6::zero)
                );
        }
    }

    return lineContactCoeffs_;
}


Foam::tensor6PairListList&
Foam::multibeamFvBlockMatrix::pointContactCoeffs()
{
    // Delete old CSR matrix
    deleteDemandDrivenData(blockRowPointersPtr_);
    deleteDemandDrivenData(blockColumnIndicesPtr_);
    deleteDemandDrivenData(blockCoeffsPtr_);
    deleteDemandDrivenData(rowPointersPtr_);
    deleteDemandDrivenData(columnIndicesPtr_);
    deleteDemandDrivenData(coeffsPtr_);
    deleteDemandDrivenData(rhsPtr_);
    deleteDemandDrivenData(solutionPtr_);

    // Set coeffs to zero
    // pointContactCoeffs_.clear();
    pointContactCoeffs_.setSize
    (
        beam_.contact().pointContacts().size()
    );

    forAll(pointContactCoeffs_, pcI)
    {
        pointContactCoeffs_[pcI] =
            tensor6PairList
            (
                2,
                tensor6Pair(tensor6::zero, tensor6::zero)
            );
    }
    
    return pointContactCoeffs_;
}


const Foam::labelList& Foam::multibeamFvBlockMatrix::blockRowPointers() const
{
    if (!blockRowPointersPtr_)
    {
        makeBlockCsrMatrix();
    }

    return *blockRowPointersPtr_;
}


const Foam::labelList& Foam::multibeamFvBlockMatrix::blockColumnIndices() const
{
    if (!blockColumnIndicesPtr_)
    {
        makeBlockCsrMatrix();
    }

    return *blockColumnIndicesPtr_;
}


const Foam::scalarField& Foam::multibeamFvBlockMatrix::blockCoeffs() const
{
    if (!blockCoeffsPtr_)
    {
        makeBlockCsrMatrix();
    }
    
    return *blockCoeffsPtr_;
}


const Foam::labelList& Foam::multibeamFvBlockMatrix::rowPointers() const
{
    if (!rowPointersPtr_)
    {
        makeCsrMatrix();
    }

    return *rowPointersPtr_;
}


const Foam::labelList& Foam::multibeamFvBlockMatrix::columnIndices() const
{
    if (!columnIndicesPtr_)
    {
        makeCsrMatrix();
    }

    return *columnIndicesPtr_;
}


const Foam::scalarField& Foam::multibeamFvBlockMatrix::coeffs() const
{
    if (!coeffsPtr_)
    {
        makeCsrMatrix();
    }
    
    return *coeffsPtr_;
}

const Foam::scalarField& Foam::multibeamFvBlockMatrix::rhs() const
{
    if (!rhsPtr_)
    {
        makeBlockCsrMatrix();
    }

    return *rhsPtr_;
}


Foam::scalarField& Foam::multibeamFvBlockMatrix::rhs()
{
    if (!rhsPtr_)
    {
        makeBlockCsrMatrix();
    }

    return *rhsPtr_;
}


const Foam::scalarField& Foam::multibeamFvBlockMatrix::solution() const
{
    if (!solutionPtr_)
    {
        makeBlockCsrMatrix();
    }

    return *solutionPtr_;
}


Foam::scalarField& Foam::multibeamFvBlockMatrix::solution()
{
    if (!solutionPtr_)
    {
        makeBlockCsrMatrix();
    }

    return *solutionPtr_;
}


void Foam::multibeamFvBlockMatrix::retriveLagrangeMultipliers()
{
    label k = nCells()*blockSize();
    for
    (
        Map<FixedList<label,3> >::const_iterator iter =
            activeConicalPulleyContacts().begin();
        iter != activeConicalPulleyContacts().end();
        ++iter
    )
    {
        label bI = (*iter)[0];
        label segI = (*iter)[1];
        label pulleyI = (*iter)[2];

        const vectorField& n =
            beam().contact().conicalPulleyContacts()[bI][segI][pulleyI]
           .normalDirection();

        vectorField& fn =
            beam_.contact().conicalPulleyContacts()[bI][segI][pulleyI]
           .normalContactForce();

        Info << "solution" << endl;
        Info << bI << ", " << segI << ", " << pulleyI << endl;
        for (label cpI=0; cpI<2; cpI++)
        {
            Info << "*" << (fn[cpI] & n[cpI]) << "*";
            fn[cpI] += n[cpI]*solution()[k];
            Info << solution()[k] << " ";
            Info << ";;" << (fn[cpI] & n[cpI]) << ";;" << endl;
            k++;
        }
        Info << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

