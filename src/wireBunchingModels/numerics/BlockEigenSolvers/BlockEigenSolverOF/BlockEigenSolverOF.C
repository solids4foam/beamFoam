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
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BlockEigenSolverOF.H"
#include "denseMatrixHelperFunctions.H"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSolverOF, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::BlockEigenSolverOF::convertFoamMatrixToEigenMatrix
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei,
    Eigen::SparseMatrix<scalar>& A
)
{
    // A monolithic rigid body adds one 6x6 block row and column after the
    // beam cells: body translation, then body rotation
    const label nRows = 6*nBlocks();

    std::vector<Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(36*(d.size() + l.size() + u.size() + 1));

    // -------------------------------------------------------------------------
    // diagonal
    // -------------------------------------------------------------------------

    label globalRowI = 0;

    forAll(d, cellI)
    {
        const scalarSquareMatrix& curD = d[cellI];

        for (label localRowI = 0; localRowI < 6; ++localRowI)
        {
            for (label localColI = 0; localColI < 6; ++localColI)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalRowI + localRowI,
                        globalRowI + localColI,
                        curD(localRowI, localColI)
                    )
                );
            }
        }

        globalRowI += 6;
    }

    // -------------------------------------------------------------------------
    // monolithic rigid body
    // -------------------------------------------------------------------------

    // The attachment cell's force and moment rows depend on the body
    // increments through the eliminated attachment face, and the body force
    // and moment rows depend on the body and on the attachment cell
    if (monolithicCoupling_.active)
    {
        const RigidBodyMonolithicCoupling& c = monolithicCoupling_;
        const label rbRow = 6*d.size();
        const label beamRow = 6*c.attachmentCell;

        // Without rotation the body's rotation rows are identity with zero
        // right-hand side
        if (!c.rotationActive)
        {
            for (label i = 0; i < 3; ++i)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>(rbRow + 3 + i, rbRow + 3 + i, 1.0)
                );
            }
        }

        for (label rowI = 0; rowI < 3; ++rowI)
        {
            for (label colI = 0; colI < 3; ++colI)
            {
                // Body translation columns of the beam rows
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        beamRow + rowI, rbRow + colI,
                        c.beamWRowCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        beamRow + 3 + rowI, rbRow + colI,
                        c.beamThetaRowCoeff(rowI, colI)
                    )
                );

                // Body translation rows
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + rowI, rbRow + colI,
                        c.bodyCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + rowI, beamRow + colI,
                        c.bodyWCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + rowI, beamRow + 3 + colI,
                        c.bodyThetaCoeff(rowI, colI)
                    )
                );

                if (!c.rotationActive)
                {
                    continue;
                }

                // Body rotation columns of the beam rows and the body
                // translation rows
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        beamRow + rowI, rbRow + 3 + colI,
                        c.beamWRowRotCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        beamRow + 3 + rowI, rbRow + 3 + colI,
                        c.beamThetaRowRotCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + rowI, rbRow + 3 + colI,
                        c.bodyTransRotCoeff(rowI, colI)
                    )
                );

                // Body rotation rows
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + 3 + rowI, rbRow + colI,
                        c.bodyRotTransCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + 3 + rowI, rbRow + 3 + colI,
                        c.bodyRotCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + 3 + rowI, beamRow + colI,
                        c.bodyRotWCoeff(rowI, colI)
                    )
                );
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        rbRow + 3 + rowI, beamRow + 3 + colI,
                        c.bodyRotThetaCoeff(rowI, colI)
                    )
                );
            }
        }
    }

    // -------------------------------------------------------------------------
    // off-diagonal
    // -------------------------------------------------------------------------

    forAll(u, faceI)
    {
        const scalarSquareMatrix& curU = u[faceI];
        const scalarSquareMatrix& curL = l[faceI];

        const label globalRow = 6*own[faceI];
        const label globalCol = 6*nei[faceI];

        // upper
        for (label localRowI = 0; localRowI < 6; ++localRowI)
        {
            for (label localColI = 0; localColI < 6; ++localColI)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalRow + localRowI,
                        globalCol + localColI,
                        curU(localRowI, localColI)
                    )
                );
            }
        }

        // lower
        for (label localRowI = 0; localRowI < 6; ++localRowI)
        {
            for (label localColI = 0; localColI < 6; ++localColI)
            {
                coefficients.push_back
                (
                    Eigen::Triplet<scalar>
                    (
                        globalCol + localRowI,
                        globalRow + localColI,
                        curL(localRowI, localColI)
                    )
                );
            }
        }
    }

    // Insert triplets into the matrix
    A.resize(nRows, nRows);
    A.setFromTriplets(coefficients.begin(), coefficients.end());
    // Compressing matrix is meant to help performance
    A.makeCompressed();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BlockEigenSolverOF::BlockEigenSolverOF
(
    const Field<scalarSquareMatrix>& d,
    const Field<scalarSquareMatrix>& l,
    const Field<scalarSquareMatrix>& u,
    const labelList& own,
    const labelList& nei
)
:
    d_(d),
    l_(l),
    u_(u),
    own_(own),
    nei_(nei),
    monolithicCoupling_(),
    bodyTranslationIncrement_(vector::zero),
    bodyRotationIncrement_(vector::zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::BlockEigenSolverOF::solve
(
    Foam::Field<Foam::scalarRectangularMatrix>& foamX,
    const Foam::Field<Foam::scalarRectangularMatrix>& foamB
)
{
    Eigen::SparseMatrix<scalar> A;
    convertFoamMatrixToEigenMatrix(d_, l_, u_, own_, nei_, A);

    const label nRows = A.rows();

    // First row of the monolithic rigid body block
    const label rigidStart = 6*d_.size();

    // -------------------------------------------------------------------------
    // Build RHS vector
    // -------------------------------------------------------------------------

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> b(nRows);
    label index = 0;

    forAll(foamB, cellI)
    {
        b(index++) = foamB[cellI](0,0);
        b(index++) = foamB[cellI](1,0);
        b(index++) = foamB[cellI](2,0);
        b(index++) = foamB[cellI](3,0);
        b(index++) = foamB[cellI](4,0);
        b(index++) = foamB[cellI](5,0);
    }

    // Monolithic rigid body: minus the body's force and moment residuals
    if (monolithicCoupling_.active)
    {
        for (label i = 0; i < 3; ++i)
        {
            b(rigidStart + i) = monolithicCoupling_.bodySource[i];
            b(rigidStart + 3 + i) =
                monolithicCoupling_.rotationActive
              ? monolithicCoupling_.bodyRotSource[i]
              : 0;
        }
    }

    // -------------------------------------------------------------------------
    // Build initial guess
    // -------------------------------------------------------------------------

    Eigen::Matrix<scalar, Eigen::Dynamic, 1> x(nRows);
    index = 0;

    forAll(foamX, cellI)
    {
        x(index++) = foamX[cellI](0,0);
        x(index++) = foamX[cellI](1,0);
        x(index++) = foamX[cellI](2,0);
        x(index++) = foamX[cellI](3,0);
        x(index++) = foamX[cellI](4,0);
        x(index++) = foamX[cellI](5,0);
    }

    // Body increments start from zero
    for (label i = rigidStart; i < nRows; ++i)
    {
        x(i) = 0.0;
    }

    // -------------------------------------------------------------------------
    // Calculate initial residual
    // -------------------------------------------------------------------------

    const label nCells = nBlocks();
    scalar initialResidual = 0.0;

    {
        Eigen::Matrix<scalar, Eigen::Dynamic, 1> Ax(nRows);
        Ax = A*x;

        Field<scalarRectangularMatrix> foamAx
        (
            nCells,
            scalarRectangularMatrix(6, 1, 0.0)
        );

        Field<scalarRectangularMatrix> foamRhs
        (
            nCells,
            scalarRectangularMatrix(6, 1, 0.0)
        );

        label k = 0;
        for (label i = 0; i < nCells; ++i)
        {
            for (label j = 0; j < 6; ++j)
            {
                foamAx[i](j,0) = Ax[k];
                foamRhs[i](j,0) = b[k];
                ++k;
            }
        }

        Field<scalarRectangularMatrix> blockR(foamRhs - foamAx);

        initialResidual = sqrt(sum(magSqr(blockR)));
    }

    // -------------------------------------------------------------------------
    // Solve
    // -------------------------------------------------------------------------

    typedef enum
    {
        SparseLU,
        BiCGSTAB,
        GMRES,
        DGMRES,
        MINRES
    } solvers;

    solvers sol = SparseLU;

    switch (sol)
    {
        case SparseLU:
        {
            Eigen::SparseLU
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::COLAMDOrdering<int>
            > solver(A);

            x = solver.solve(b);
            break;
        }

        case BiCGSTAB:
        {
            Eigen::BiCGSTAB
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::IncompleteLUT<scalar>
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        case GMRES:
        {
            Eigen::GMRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::IncompleteLUT<scalar>
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        case DGMRES:
        {
            Eigen::DGMRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::IncompleteLUT<scalar>
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        case MINRES:
        {
            Eigen::MINRES
            <
                Eigen::SparseMatrix<scalar>,
                Eigen::Lower|Eigen::Upper,
                Eigen::IdentityPreconditioner
            > solver;

            solver.compute(A);
            x = solver.solve(b);
            break;
        }

        default:
            FatalErrorIn
            (
                "Foam::scalar Foam::BlockEigenSolverOF::solve(...)"
            )   << "Undefined Eigen solver."
                << abort(FatalError);
    }

    // -------------------------------------------------------------------------
    // Copy solved beam unknowns back to foamX
    // -------------------------------------------------------------------------

    index = 0;
    forAll(foamX, cellI)
    {
        foamX[cellI](0,0) = x(index++);
        foamX[cellI](1,0) = x(index++);
        foamX[cellI](2,0) = x(index++);
        foamX[cellI](3,0) = x(index++);
        foamX[cellI](4,0) = x(index++);
        foamX[cellI](5,0) = x(index++);
    }

    // Monolithic rigid body: Newton increments of the body displacement and
    // of its step rotation vector
    if (monolithicCoupling_.active)
    {
        for (label i = 0; i < 3; ++i)
        {
            bodyTranslationIncrement_[i] = x(rigidStart + i);
            bodyRotationIncrement_[i] = x(rigidStart + 3 + i);
        }
    }
    else
    {
        bodyTranslationIncrement_ = vector::zero;
        bodyRotationIncrement_ = vector::zero;
    }

    return initialResidual;
}


// ************************************************************************* //
