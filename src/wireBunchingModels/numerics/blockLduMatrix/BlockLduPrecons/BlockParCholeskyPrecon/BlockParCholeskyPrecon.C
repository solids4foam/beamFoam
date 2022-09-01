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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "BlockParCholeskyPrecon.H"
#include "BlockLduInterfaceFieldPtrsList.H"
#include "ProcessorBlockLduInterfaceField.H"
#include "processorLduInterfaceField.H"
#include "processorFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockParCholeskyPrecon<Type>::calcPreconDiag()
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    // Precondition the diagonal

    if (this->matrix_.symmetric())
    {
        // Get interface list
        const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces =
            this->matrix_.interfaces();

        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asScalar(),
                            this->matrix_.coupleUpper()[patchI].asScalar(),
                            this->matrix_.coupleLower()[patchI].asScalar()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asLinear(),
                            this->matrix_.coupleUpper()[patchI].asLinear(),
                            this->matrix_.coupleLower()[patchI].asLinear()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asSquare(),
                            this->matrix_.coupleLower()[patchI].asSquare()
                        );
                    }
                }

                // Do core matrix
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asLinear(),
                            this->matrix_.coupleUpper()[patchI].asScalar(),
                            this->matrix_.coupleLower()[patchI].asScalar()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asLinear(),
                            this->matrix_.coupleUpper()[patchI].asLinear(),
                            this->matrix_.coupleLower()[patchI].asLinear()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asSquare(),
                            this->matrix_.coupleLower()[patchI].asSquare()
                        );
                    }
                }

                // Do core matrix
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asScalar(),
                            this->matrix_.coupleLower()[patchI].asScalar()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asLinear(),
                            this->matrix_.coupleLower()[patchI].asLinear()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asSquare(),
                            this->matrix_.coupleLower()[patchI].asSquare()
                        );
                    }
                }

                // Do core matrix
                // Transpose multiplication
                diagMultiplyCoeffT
                (
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
    }
    else // Asymmetric matrix
    {
        // Get interface list
        const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces =
            this->matrix_.interfaces();

        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asScalar(),
                            this->matrix_.coupleUpper()[patchI].asScalar(),
                            this->matrix_.coupleLower()[patchI].asScalar()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asLinear(),
                            this->matrix_.coupleUpper()[patchI].asLinear(),
                            this->matrix_.coupleLower()[patchI].asLinear()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asSquare(),
                            this->matrix_.coupleLower()[patchI].asSquare()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asLinear(),
                            this->matrix_.coupleUpper()[patchI].asScalar(),
                            this->matrix_.coupleLower()[patchI].asScalar()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asLinear(),
                            this->matrix_.coupleUpper()[patchI].asLinear(),
                            this->matrix_.coupleLower()[patchI].asLinear()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asSquare(),
                            this->matrix_.coupleLower()[patchI].asSquare()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare()
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asScalar(),
                            this->matrix_.coupleLower()[patchI].asScalar()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                // Do coupled interfaces
                forAll (interfaces, patchI)
                {
                    if (interfaces.set(patchI))
                    {
                        // Get face-cells addressing
                        const unallocLabelList& fc =
                            interfaces[patchI].coupledInterface().faceCells();

                        diagInterfaceMultiply
                        (
                            fc,
                            preconDiag_.asSquare(),
                            this->matrix_.coupleUpper()[patchI].asLinear(),
                            this->matrix_.coupleLower()[patchI].asLinear()
                        );
                    }
                }

                // Do core matrix
                diagMultiply
                (
                    preconDiag_.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear()
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Pout << mag(preconDiag_.asSquare()) << endl;

                
                typedef typename BlockCoeff<Type>::squareType coeffType;
                
                // Get addressing
                const unallocLabelList& upperAddr =
                    this->matrix_.lduAddr().upperAddr();
                const unallocLabelList& lowerAddr =
                    this->matrix_.lduAddr().lowerAddr();

                lduInterfaceFieldPtrsList interfaces_
                (
                    interfaces.size()
                );
                forAll(interfaces, intI)
                {
                    if (interfaces.set(intI))
                    {
                        interfaces_.set(intI, interfaces(intI));
                    }
                }

                // Get list of coefficients for internal equations
                // (faces where at least
                // one cell is an internal cell, i.e. not a boundary cell)
                const unallocLabelList& internalCoeffs =
                    this->matrix_.lduAddr().internalEqnCoeffs
                    (
                        interfaces_
                    );

                // Get list of coefficients for internal equations
                // where the upper/lower
                // need to be flipped (faces where a face cell is the owner of a
                // neighbouring internal cell)
                const unallocLabelList& flippedInternalCoeffs =
                    this->matrix_.lduAddr().flippedInternalEqnCoeffs
                    (
                        interfaces_
                    );

                // Create multiplication function object
                typename BlockCoeff<Type>::multiply mult;
    
                // STAGE 1: Perform factorisation for internal equations only

                // Loop through coeffs (internal faces) for internal equations
                // (cells) that do not need a flip
                forAll (internalCoeffs, icI)
                {
                    // Get coefficient label
                    const label& coeffI = internalCoeffs[icI];

                    // Update upper (neighbour of the face)
                    preconDiag_.asSquare()[upperAddr[coeffI]] -=
                        mult.tripleProduct
                        (
                            LowerCoeff.asSquare()[coeffI],
                            preconDiag_.asSquare()[lowerAddr[coeffI]],
                            UpperCoeff.asSquare()[coeffI]
                        );
                }

                // Loop through remaining coeffs (internal faces)
                // for internal equations
                // (cells) that need a flip
                forAll (flippedInternalCoeffs, ficI)
                {
                    // Get coefficient label
                    const label& coeffI = flippedInternalCoeffs[ficI];

                    // Update lower (owner of the face)
                    preconDiag_.asSquare()[lowerAddr[coeffI]] -=
                        mult.tripleProduct
                        (
                            UpperCoeff.asSquare()[coeffI],
                            preconDiag_.asSquare()[upperAddr[coeffI]],
                            LowerCoeff.asSquare()[coeffI]
                        );

                    Pout << coeffI << endl;
                }


                // STAGE 2: Factorisation across coupled boundaries

                // Coupled interface update: note, ordering is important
                if (true)
                {
                    // Placeholder field for inverted diagonal at face cells that needs
                    // to be populated just before sending the data
                    Field<coeffType> diag(preconDiag_.size(), pTraits<coeffType>::zero);

                    forAll (interfaces, intI)
                    {
                        if (interfaces.set(intI))
                        {
                            // Get boundary coefficients for this interface (faces
                            // between two boundary cells for this patch)
                            const dynamicLabelList& boundaryCoeffs =
                                this->matrix_.lduAddr().boundaryEqnCoeffs(interfaces_, intI);

                            if (interfaces_[intI].coupledInterface().master())
                            {
                                // This is master: need to take into account boundary
                                // coeff contributions first
                                forAll (boundaryCoeffs, bcI)
                                {
                                    // Get coefficient label
                                    const label& coeffI = boundaryCoeffs[bcI];

                                    preconDiag_.asSquare()[upperAddr[coeffI]] -=
                                        mult.tripleProduct
                                        (
                                            LowerCoeff.asSquare()[coeffI],
                                            preconDiag_.asSquare()[lowerAddr[coeffI]],
                                            UpperCoeff.asSquare()[coeffI]
                                        );
                                }

                                // Calculate inverse diagonal for face cells of this
                                // boundary since master needs to send them
                                const unallocLabelList& fc =
                                    interfaces[intI].coupledInterface().faceCells();

                                forAll (fc, fcI)
                                {
                                    const label& cellI = fc[fcI];

                                    diag[cellI] = preconDiag_.asSquare()[cellI];
                                        // 1.0/preconDiag_[cellI];
                                }

                                // Perform coupled boundary update for this interface

                                if (isA<processorLduInterfaceField>(interfaces[intI]))
                                {
                                    const processorLduInterface& procPatch =
                                        refCast<const processorLduInterface>
                                        (
                                            this->matrix_.mesh().interfaces()[intI]
                                        );

                                    Field<coeffType> sendDiag(diag, fc);

                                    // Send inverted diagonal
                                    procPatch.send
                                    (
                                        Pstream::blocking,
                                        sendDiag
                                    );
                                }
                                else
                                {
                                    FatalErrorIn
                                    (
                                        "void Foam::BlockParCholeskyPrecon<Type>::calcPreconDiag()"
                                    )   << "Parallel ILU0 precondicioner not implemented for "
                                        << interfaces[intI].type() << exit(FatalIOError);
                                }
                            }
                            else
                            {
                                const unallocLabelList& fc =
                                    interfaces[intI].coupledInterface().faceCells();
                                
                                // This is slave: need to perform coupled boundary
                                // update for this interface first

                                if (isA<processorLduInterfaceField>(interfaces[intI]))
                                {
                                    const processorLduInterface& procPatch =
                                        refCast<const processorLduInterface>
                                        (
                                            this->matrix_.mesh().interfaces()[intI]
                                        );
                                    
                                    Field<coeffType> receiveDiag
                                    (
                                        fc.size(),
                                        pTraits<coeffType>::zero
                                    );
                                    
                                    procPatch.receive<coeffType>
                                    (
                                        Pstream::blocking,
                                        receiveDiag
                                    );

                                    forAll(fc, fI)
                                    {
                                        preconDiag_.asSquare()[fc[fI]] -=
                                            mult.tripleProduct
                                            (
                                                this->matrix_.coupleUpper()[intI].asSquare()[fI],
                                                receiveDiag[fI],
                                                this->matrix_.coupleLower()[intI].asSquare()[fI]
                                            );
                                    }
                                }
                                else
                                {
                                    FatalErrorIn
                                    (
                                        "void Foam::BlockParCholeskyPrecon<Type>::calcPreconDiag()"
                                    )   << "Parallel ILU0 precondicioner not implemented for "
                                        << interfaces[intI].type() << exit(FatalIOError);
                                }

                                // Update boundary coeffs after downstream coupled
                                // boundary contribution is taken into account
                                forAll (boundaryCoeffs, bcI)
                                {
                                    // Get coefficient label
                                    const label& coeffI = boundaryCoeffs[bcI];

                                    preconDiag_.asSquare()[upperAddr[coeffI]] -=
                                        mult.tripleProduct
                                        (
                                            LowerCoeff.asSquare()[coeffI],
                                            preconDiag_.asSquare()[lowerAddr[coeffI]],
                                            UpperCoeff.asSquare()[coeffI]
                                        );
                                }
                            }
                        }
                    }
                }
                
                // // Do coupled interfaces
                // forAll (interfaces, patchI)
                // {
                //     if (interfaces.set(patchI))
                //     {
                //         // Get face-cells addressing
                //         const unallocLabelList& fc =
                //             interfaces[patchI].coupledInterface().faceCells();

                //         diagInterfaceMultiply
                //         (
                //             fc,
                //             preconDiag_.asSquare(),
                //             this->matrix_.coupleUpper()[patchI].asSquare(),
                //             this->matrix_.coupleLower()[patchI].asSquare()
                //         );
                //     }
                // }

                // // Do core matrix
                // diagMultiply
                // (
                //     preconDiag_.asSquare(),
                //     LowerCoeff.asSquare(),
                //     UpperCoeff.asSquare()
                // );
            }
        }
    }

    // Invert the diagonal
    preconDiag_ = inv(preconDiag_);

     // Pout << mag(preconDiag_.asSquare()) << endl;
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::diagMultiply
(
    Field<DiagType>& dDiag,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                upper[coeffI],
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::diagMultiplyCoeffT
(
    Field<DiagType>& dDiag,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                upper[coeffI].T(),        // Upper coefficient transposed
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::diagMultiply
(
    Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper
)
{
    // Precondition the diagonal

    // Get addressing
    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (upper, coeffI)
    {
        dDiag[upperAddr[coeffI]] -=
            mult.tripleProduct
            (
                lower[coeffI],
                dDiag[lowerAddr[coeffI]],
                upper[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::diagInterfaceMultiply
(
    const unallocLabelList& fc,
    Field<DiagType>& dDiag,
    const Field<ULType>& bouCoeffs,
    const Field<ULType>& intCoeffs
)
{
    // Precondition the diagonal for the coupled interface

    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (fc, coeffI)
    {
        // Note: possible sign issue.  HJ and VV, 19/Jun/2017
        dDiag[fc[coeffI]] +=
            mult.tripleProduct
            (
                intCoeffs[coeffI],
                dDiag[fc[coeffI]],
                bouCoeffs[coeffI]
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::ILUmultiply
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    forAll (upper, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI], x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::ILUmultiplyCoeffT
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();

    forAll (upper, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI].T(), x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::ILUmultiply
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = this->matrix_.lduAddr().losortAddr();

    label losortCoeff;

    forAll (lower, coeffI)
    {
        losortCoeff = losortAddr[coeffI];

        x[upperAddr[losortCoeff]] -=
            mult
            (
                dDiag[upperAddr[losortCoeff]],
                mult(lower[losortCoeff], x[lowerAddr[losortCoeff]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        x[lowerAddr[coeffI]] -=
            mult
            (
                dDiag[lowerAddr[coeffI]],
                mult(upper[coeffI], x[upperAddr[coeffI]])
            );
    }
}


template<class Type>
template<class DiagType, class ULType>
void Foam::BlockParCholeskyPrecon<Type>::ILUmultiplyTranspose
(
    Field<Type>& x,
    const Field<DiagType>& dDiag,
    const Field<ULType>& lower,
    const Field<ULType>& upper,
    const Field<Type>& b
) const
{
    // Create multiplication function object
    typename BlockCoeff<Type>::multiply mult;

    forAll (x, i)
    {
        x[i] = mult(dDiag[i], b[i]);
    }

    const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
    const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();
    const unallocLabelList& losortAddr = this->matrix_.lduAddr().losortAddr();

    label losortCoeff;

    //HJ Not sure if the coefficient itself needs to be transposed.
    // HJ, 30/Oct/2007
    forAll (lower, coeffI)
    {
        x[upperAddr[coeffI]] -=
            mult
            (
                dDiag[upperAddr[coeffI]],
                mult(upper[coeffI], x[lowerAddr[coeffI]])
            );
    }

    forAllReverse (upper, coeffI)
    {
        losortCoeff = losortAddr[coeffI];

        x[lowerAddr[losortCoeff]] -=
            mult
            (
                dDiag[lowerAddr[losortCoeff]],
                mult(lower[losortCoeff], x[upperAddr[losortCoeff]])
            );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockParCholeskyPrecon<Type>::BlockParCholeskyPrecon
(
    const BlockLduMatrix<Type>& matrix
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag())
{
    this->calcPreconDiag();
}


template<class Type>
Foam::BlockParCholeskyPrecon<Type>::BlockParCholeskyPrecon
(
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduPrecon<Type>(matrix),
    preconDiag_(matrix.diag())
{
    calcPreconDiag();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockParCholeskyPrecon<Type>::precondition
(
    Field<Type>& x,
    const Field<Type>& b
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asScalar(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asLinear(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // Multiplication using transposed upper square coefficient
                ILUmultiplyCoeffT
                (
                    x,
                    preconDiag_.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asScalar(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    b
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiply
                (
                    x,
                    preconDiag_.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    b
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                // typedef typename BlockCoeff<Type>::squareType coeffType;
                
                // Create multiplication function object
                typename BlockCoeff<Type>::multiply mult;

                forAll (x, i)
                {
                    x[i] = mult(preconDiag_.asSquare()[i], b[i]);
                }

                const unallocLabelList& upperAddr = this->matrix_.lduAddr().upperAddr();
                const unallocLabelList& lowerAddr = this->matrix_.lduAddr().lowerAddr();
                const unallocLabelList& losortAddr = this->matrix_.lduAddr().losortAddr();

                // Get interface list
                const typename BlockLduInterfaceFieldPtrsList<Type>::Type& interfaces =
                    this->matrix_.interfaces();

                const TypeCoeffField& LowerCoeff = this->matrix_.lower();
                const TypeCoeffField& UpperCoeff = this->matrix_.upper();

                lduInterfaceFieldPtrsList interfaces_
                (
                    interfaces.size()
                );
                forAll(interfaces, intI)
                {
                    if (interfaces.set(intI))
                    {
                        interfaces_.set(intI, interfaces(intI));
                    }
                }

                // Get list of coefficients for internal equations
                // (faces where at least
                // one cell is an internal cell, i.e. not a boundary cell)
                const unallocLabelList& internalCoeffs =
                    this->matrix_.lduAddr().internalEqnCoeffs
                    (
                        interfaces_
                    );

                // Get list of coefficients for internal equations
                // where the upper/lower
                // need to be flipped (faces where a face cell is the owner of a
                // neighbouring internal cell)
                const unallocLabelList& flippedInternalCoeffs =
                    this->matrix_.lduAddr().flippedInternalEqnCoeffs
                    (
                        interfaces_
                    );

                // STAGE 1: Forward substitution

                // Register helper variables
                register label coeffI = -1, losortCoeffI = -1, u = -1, l = -1;

                // Loop through coeffs (internal faces) for internal equations (cells)
                // that do not need a flip
                forAll (internalCoeffs, icI)
                {
                    // Get index for this coefficient
                    losortCoeffI = losortAddr[internalCoeffs[icI]];

                    // Get lower/upper
                    l = lowerAddr[losortCoeffI];
                    u = upperAddr[losortCoeffI];

                    // Update upper (neighbour of the face)
                    x[u] -= //preconDiag_[u]*lower[losortCoeffI]*x[l];
                        mult
                        (
                            preconDiag_.asSquare()[u],
                            mult(LowerCoeff.asSquare()[losortCoeffI], x[l])
                        );
                }

                // Loop through remaining coeffs (internal faces) for internal equations
                // (cells) that need a flip
                forAll (flippedInternalCoeffs, ficI)
                {
                    // Get index for this coefficient
                    losortCoeffI = losortAddr[flippedInternalCoeffs[ficI]];

                    // Get lower/upper
                    l = lowerAddr[losortCoeffI];
                    u = upperAddr[losortCoeffI];
                    
                    // Update lower (owner of the face)
                    x[l] -= //preconDiag_[l]*upper[losortCoeffI]*x[u];
                        mult
                        (
                            preconDiag_.asSquare()[l],
                            mult(UpperCoeff.asSquare()[losortCoeffI], x[u])
                        );
                }

                
                // Coupled interface update for the forward solve

                // Note: ordering is important
                forAll (interfaces, intI)
                {
                    if (interfaces.set(intI))
                    {
                        // Get boundary coefficients for this interface (faces between
                        // two boundary cells for this patch)
                        const dynamicLabelList& boundaryCoeffs =
                            this->matrix_.lduAddr().boundaryEqnCoeffs(interfaces_, intI);

                        if (interfaces[intI].coupledInterface().master())
                        {
                            // This is master: need to take into account boundary coeff
                            // contributions first
                            
                            forAll (boundaryCoeffs, bcI)
                            {
                                // Get index for this coefficient
                                losortCoeffI = losortAddr[boundaryCoeffs[bcI]];

                                // Get lower/upper
                                l = lowerAddr[losortCoeffI];
                                u = upperAddr[losortCoeffI];
                                
                                x[u] -= //preconDiag_[u]*lower[coeffI]*x[l];
                                    mult
                                    (
                                        preconDiag_.asSquare()[u],
                                        mult(LowerCoeff.asSquare()[losortCoeffI], x[l])
                                    );
                            }

                            const unallocLabelList& fc =
                                interfaces[intI].coupledInterface().faceCells();
                                
                            if (isA<processorLduInterfaceField>(interfaces[intI]))
                            {
                                const processorLduInterface& procPatch =
                                    refCast<const processorLduInterface>
                                    (
                                        this->matrix_.mesh().interfaces()[intI]
                                    );

                                Field<Type> sendx(x, fc);

                                // Send x
                                procPatch.send
                                (
                                    Pstream::blocking,
                                    sendx
                                );
                            }
                            else
                            {
                                FatalErrorIn
                                (
                                    "void Foam::BlockParCholeskyPrecon<Type>::preconditions\n"
                                    "(\n"
                                    "    Field<Type>& x,\n"
                                    "    const Field<Type>& b\n"
                                    ") const\n"
                                )   << "Parallel ILU0 precondicioner not implemented for "
                                    << interfaces[intI].type() << exit(FatalIOError);
                            }
                        }
                        else
                        {
                            // This is slave: need to perform coupled boundary update
                            // for this interface first

                            const unallocLabelList& fc =
                                interfaces[intI].coupledInterface().faceCells();
                            
                            if (isA<processorLduInterfaceField>(interfaces[intI]))
                            {
                                const processorLduInterface& procPatch =
                                    refCast<const processorLduInterface>
                                    (
                                        this->matrix_.mesh().interfaces()[intI]
                                    );
                                    
                                Field<Type> receivex
                                (
                                    fc.size(),
                                    pTraits<Type>::zero
                                );
                                    
                                procPatch.receive<Type>
                                (
                                    Pstream::blocking,
                                    receivex
                                );

                                forAll(fc, fI)
                                {
                                    x[fc[fI]] += // coupleUpper has got changed sign 
                                        mult
                                        (
                                            preconDiag_.asSquare()[fc[fI]],
                                            mult
                                            (
                                                this->matrix_.coupleUpper()[intI].asSquare()[fI],
                                                receivex[fI]
                                            )
                                        );
                                }
                            }
                            else
                            {
                                FatalErrorIn
                                (
                                    "void Foam::BlockParCholeskyPrecon<Type>::preconditions\n"
                                    "(\n"
                                    "    Field<Type>& x,\n"
                                    "    const Field<Type>& b\n"
                                    ") const\n"
                                )   << "Parallel ILU0 precondicioner not implemented for "
                                    << interfaces[intI].type() << exit(FatalIOError);
                            }

                            // Update boundary coeffs after downstream coupled boundary
                            // contribution is taken into account
                            forAll (boundaryCoeffs, bcI)
                            {
                                // Get index for this coefficient
                                losortCoeffI = losortAddr[boundaryCoeffs[bcI]];

                                // Get lower/upper
                                l = lowerAddr[losortCoeffI];
                                u = upperAddr[losortCoeffI];

                                x[u] -= //preconDiag_[u]*lower[coeffI]*x[l];
                                    mult
                                    (
                                        preconDiag_.asSquare()[u],
                                        mult(LowerCoeff.asSquare()[losortCoeffI], x[l])
                                    );
                            }
                        }
                    }
                }

                // STAGE 2: Backward substitution

                // Note: ordering is important and since this is a backward substitution
                // part, master/slave have exchanged roles and flipped order of
                // execution of coupled boundary update and boundaryCoeffs update
                forAll (interfaces, intI)
                {
                    if (interfaces.set(intI))
                    {
                        // Get boundary coefficients for this interface (faces between
                        // two boundary cells for this patch)
                        const dynamicLabelList& boundaryCoeffs =
                            this->matrix_.lduAddr().boundaryEqnCoeffs(interfaces_, intI);

                        if (interfaces[intI].coupledInterface().master())
                        {
                            // This is master: need to update coupled boundaries first.
                            // Note: in the backward substitution, the information
                            // propagates from slave to master (not from master to
                            // slave). Hence, coupled boundary update needs to come
                            // first on the master side!
                            
                            const unallocLabelList& fc =
                                interfaces[intI].coupledInterface().faceCells();
                            
                            if (isA<processorLduInterfaceField>(interfaces[intI]))
                            {
                                const processorLduInterface& procPatch =
                                    refCast<const processorLduInterface>
                                    (
                                        this->matrix_.mesh().interfaces()[intI]
                                    );
                                    
                                Field<Type> receivex
                                (
                                    fc.size(),
                                    pTraits<Type>::zero
                                );
                                    
                                procPatch.receive<Type>
                                (
                                    Pstream::blocking,
                                    receivex
                                );

                                forAll(fc, fI)
                                {
                                    x[fc[fI]] += // coupleUpper has got changed sign 
                                        mult
                                        (
                                            preconDiag_.asSquare()[fc[fI]],
                                            mult
                                            (
                                                this->matrix_.coupleUpper()[intI].asSquare()[fI],
                                                receivex[fI]
                                            )
                                        );
                                }
                            }
                            else
                            {
                                FatalErrorIn
                                (
                                    "void Foam::BlockParCholeskyPrecon<Type>::preconditions\n"
                                    "(\n"
                                    "    Field<Type>& x,\n"
                                    "    const Field<Type>& b\n"
                                    ") const\n"
                                )   << "Parallel ILU0 precondicioner not implemented for "
                                    << interfaces[intI].type() << exit(FatalIOError);
                            }
                            
                            // Update boundary coeffs after upstream coupled boundary
                            // contribution is taken into account. Note: reverse loop
                            forAllReverse(boundaryCoeffs, bcI)
                            {
                                // Get coefficient label
                                coeffI = boundaryCoeffs[bcI];

                                // Get lower/upper
                                l = lowerAddr[coeffI];
                                u = upperAddr[coeffI];

                                x[l] -= //preconDiag_[l]*upper[coeffI]*x[u];
                                    mult
                                    (
                                        preconDiag_.asSquare()[l],
                                        mult(UpperCoeff.asSquare()[coeffI], x[u])
                                    );
                            }
                        }
                        else
                        {
                            // This is slave: need to take into account boundary coeff
                            // contributions first

                            forAllReverse (boundaryCoeffs, bcI)
                            {
                                // Get coefficient label
                                coeffI = boundaryCoeffs[bcI];

                                // Get lower/upper
                                l = lowerAddr[coeffI];
                                u = upperAddr[coeffI];

                                x[l] -= //preconDiag_[l]*upper[coeffI]*x[u];
                                    mult
                                    (
                                        preconDiag_.asSquare()[l],
                                        mult(UpperCoeff.asSquare()[coeffI], x[u])
                                    );
                            }

                            const unallocLabelList& fc =
                                interfaces[intI].coupledInterface().faceCells();

                            if (isA<processorLduInterfaceField>(interfaces[intI]))
                            {
                                const processorLduInterface& procPatch =
                                    refCast<const processorLduInterface>
                                    (
                                        this->matrix_.mesh().interfaces()[intI]
                                    );

                                Field<Type> sendx(x, fc);

                                // Send x
                                procPatch.send
                                (
                                    Pstream::blocking,
                                    sendx
                                );
                            }
                            else
                            {
                                FatalErrorIn
                                (
                                    "void Foam::BlockParCholeskyPrecon<Type>::preconditions\n"
                                    "(\n"
                                    "    Field<Type>& x,\n"
                                    "    const Field<Type>& b\n"
                                    ") const\n"
                                )   << "Parallel ILU0 precondicioner not implemented for "
                                    << interfaces[intI].type() << exit(FatalIOError);
                            }
                        }
                    }
                }
                
                // Loop through coeffs (internal faces) for internal equations (cells)
                // that do not need a flip, after coupled boundary update
                forAllReverse (internalCoeffs, icI)
                {
                    // Get index for this coefficient
                    coeffI = internalCoeffs[icI];

                    // Get lower/upper
                    l = lowerAddr[coeffI];
                    u = upperAddr[coeffI];

                    // Update lower (owner of the face)
                    x[l] -= //preconDiag_[l]*upper[coeffI]*x[u];
                        mult
                        (
                            preconDiag_.asSquare()[l],
                            mult(UpperCoeff.asSquare()[coeffI], x[u])
                        );
                }

                // Loop through remaining coeffs (internal faces) for internal equations
                // (cells) that need a flip, after coupled boundary update
                forAllReverse (flippedInternalCoeffs, ficI)
                {
                    // Get index for this coefficient
                    coeffI = flippedInternalCoeffs[ficI];

                    // Get lower/upper
                    l = lowerAddr[coeffI];
                    u = upperAddr[coeffI];

                    // Update upper (neighbour of the face)
                    x[u] -= //preconDiag_[u]*lower[coeffI]*x[l];
                        mult
                        (
                            preconDiag_.asSquare()[u],
                            mult(LowerCoeff.asSquare()[coeffI], x[u])
                        );
                }
        
                // ILUmultiply
                // (
                //     x,
                //     preconDiag_.asSquare(),
                //     LowerCoeff.asSquare(),
                //     UpperCoeff.asSquare(),
                //     b
                // );
            }
        }
    }// end axi-symmetric
}

template<class Type>
void Foam::BlockParCholeskyPrecon<Type>::preconditionT
(
    Field<Type>& xT,
    const Field<Type>& bT
) const
{
    typedef CoeffField<Type> TypeCoeffField;

    // Note: Assuming lower and upper triangle have the same active type

    if (this->matrix_.symmetric())
    {
        precondition(xT, bT);
    }
    else // Asymmetric matrix
    {
        const TypeCoeffField& LowerCoeff = this->matrix_.lower();
        const TypeCoeffField& UpperCoeff = this->matrix_.upper();

        if (preconDiag_.activeType() == blockCoeffBase::SCALAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asScalar(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::LINEAR)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asLinear(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
        else if (preconDiag_.activeType() == blockCoeffBase::SQUARE)
        {
            if (UpperCoeff.activeType() == blockCoeffBase::SCALAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asSquare(),
                    LowerCoeff.asScalar(),
                    UpperCoeff.asScalar(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::LINEAR)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asSquare(),
                    LowerCoeff.asLinear(),
                    UpperCoeff.asLinear(),
                    bT
                );
            }
            else if (UpperCoeff.activeType() == blockCoeffBase::SQUARE)
            {
                ILUmultiplyTranspose
                (
                    xT,
                    preconDiag_.asSquare(),
                    LowerCoeff.asSquare(),
                    UpperCoeff.asSquare(),
                    bT
                );
            }
        }
    }
}


template<class Type>
void Foam::BlockParCholeskyPrecon<Type>::initMatrix()
{
    preconDiag_ = this->matrix_.diag();

    this->calcPreconDiag();
}


// ************************************************************************* //
