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

Description
    Preconditioned Bi-Conjugate Gradient stabilised solver.

\*---------------------------------------------------------------------------*/

#include "myBlockBiCGStabSolver.H"
#include "addToRunTimeSelectionTable.H"
// #include<unistd.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myBlockBiCGStabSolver, 0);
    addToRunTimeSelectionTable
    (
        blockBiCGStabVector6Solver, myBlockBiCGStabSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockBiCGStabVector6Solver, myBlockBiCGStabSolver, asymMatrix
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::myBlockBiCGStabSolver::myBlockBiCGStabSolver
(
    const Foam::word& fieldName,
    const Foam::BlockLduMatrix<vector6>& matrix,
    const Foam::dictionary& dict
)
:
    blockBiCGStabVector6Solver
    (
        fieldName,
        matrix,
        dict
        ),
    preconPtr_
    (
        BlockLduPrecon<vector6>::New
        (
            matrix,
            this->dict()
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

typename Foam::BlockSolverPerformance<Foam::vector6>
Foam::myBlockBiCGStabSolver::solve
(
    Field<vector6>& x,
    const Field<vector6>& b
)
{
    // Create local references to avoid the spread this-> ugliness
    const BlockLduMatrix<vector6>& matrix = this->matrix_;

    // Prepare solver performance
    BlockSolverPerformance<vector6> solverPerf
    (
        typeName,
        this->fieldName()
    );

    vector6 norm = this->normFactor(x, b);

    // Multiplication helper
    typename BlockCoeff<vector6>::multiply mult;

    Field<vector6> p(x.size());

    // Calculate initial residual
    matrix.Amul(p, x);
    Field<vector6> r(b - p);

    solverPerf.initialResidual() = cmptDivide(gSum(cmptMag(r)),norm);
    // solverPerf.initialResidual() = gSum(cmptMag(r))/norm;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged

    if (!this->stop(solverPerf))
    {
        scalar rho = this->great_;
        scalar rhoOld = rho;

        scalar alpha = 0;
        scalar omega = this->great_;
        scalar beta;

        p = pTraits<vector6>::zero;
        Field<vector6> ph(x.size(), pTraits<vector6>::zero);
        Field<vector6> v(x.size(), pTraits<vector6>::zero);
        Field<vector6> s(x.size(), pTraits<vector6>::zero);
        Field<vector6> sh(x.size(), pTraits<vector6>::zero);
        Field<vector6> t(x.size(), pTraits<vector6>::zero);

        // Calculate transpose residual
        Field<vector6> rw(r);

        do
        {
            rhoOld = rho;

            // Update search directions
            rho = gSumProd(rw, r);

            beta = rho/rhoOld*(alpha/omega);

            // Restart if breakdown occurs
            if (rho == 0)
            {
                rw = r;
                rho = gSumProd(rw, r);

                alpha = 0;
                omega = 0;
                beta = 0;
            }

            forAll (p, i)
            {
                p[i] = r[i] + beta*p[i] - beta*omega*v[i];
            }

            preconPtr_->precondition(ph, p);

            matrix.Amul(v, ph);
            
            alpha = rho/gSumProd(rw, v);

            forAll (s, i)
            {
                s[i] = r[i] - alpha*v[i];
            }

            // Bug fix, Alexander Monakov, 11/Jul/2012
            preconPtr_->precondition(sh, s);
            matrix.Amul(t, sh);
            omega = gSumProd(t, s)/gSumProd(t, t);

            forAll (x, i)
            {
                x[i] = x[i] + alpha*ph[i] + omega*sh[i];
            }

            forAll (r, i)
            {
                r[i] = s[i] - omega*t[i];
            }

            solverPerf.finalResidual() = cmptDivide(gSum(cmptMag(r)), norm);
            // solverPerf.finalResidual() = gSum(cmptMag(r))/norm;
            solverPerf.nIterations()++;
        } while (!this->stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
