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

#include "BlockTDMSolver.H"
#include "addToRunTimeSelectionTable.H"
// #include<unistd.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockTDMSolver, 0);
    addToRunTimeSelectionTable
    (
        blockBiCGStabVector6Solver, BlockTDMSolver, symMatrix
    );

    addToRunTimeSelectionTable
    (
        blockBiCGStabVector6Solver, BlockTDMSolver, asymMatrix
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::BlockTDMSolver::BlockTDMSolver
(
    const Foam::word& fieldName,
    const Foam::BlockLduMatrix<vector6>& matrix,
    const Foam::dictionary& dict
)
:
    BlockIterativeSolver<vector6>
    (
        fieldName,
        matrix,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

typename Foam::BlockSolverPerformance<Foam::vector6>
Foam::BlockTDMSolver::solve
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
    
    vector6 norm = vector6::one; //this->normFactor(x, b);
    
    // Calculate initial residual
    {
        Field<vector6> p(x.size());
        matrix.Amul(p, x);
        Field<vector6> r(b - p);

        solverPerf.initialResidual() = cmptDivide(gSum(cmptMag(r)),norm);
        // solverPerf.finalResidual() = solverPerf.initialResidual();
    }

    // Solve
    const Field<tensor6>& l = matrix.lower().asSquare();
    const Field<tensor6>& u = matrix.upper().asSquare();
    const Field<tensor6>& d = matrix.diag().asSquare();

    // Multiplication helper
    typename BlockCoeff<vector6>::multiply mult;

    Field<tensor6> uPrime(x.size(), pTraits<tensor6>::zero);
    Field<vector6> bPrime(x.size(), pTraits<vector6>::zero);

    uPrime[0] = mult.activeTypeMultiply(inv(d[0]), u[0]);
    bPrime[0] = mult(inv(d[0]), b[0]);
    
    for (label i=1; i<u.size(); i++)
    {
        uPrime[i] =
            mult.activeTypeMultiply
            (
                inv
                (
                    d[i]
                  - mult.activeTypeMultiply(l[i-1], uPrime[i-1])
                ),
                u[i]
            );
    }

    for (label i=1; i<b.size(); i++)
    {
        bPrime[i] =
            mult
            (
                inv
                (
                    d[i]
                  - mult.activeTypeMultiply(l[i-1], uPrime[i-1])
                ),
                b[i]
              - mult(l[i-1], bPrime[i-1])
            );
    }
    
    x[x.size()-1] = bPrime[x.size()-1];

    for (label i=x.size()-2; i>=0; i--)
    {
        x[i] = bPrime[i] - mult(uPrime[i], x[i+1]);
    }
    
    // Check convergence, solve if not converged

    // if (!this->stop(solverPerf))
    // {
    //     // scalar rho = this->great_;
    //     // scalar rhoOld = rho;

    //     // scalar alpha = 0;
    //     // scalar omega = this->great_;
    //     // scalar beta;

    //     Field<vector6> uPrime(x.size(), pTraits<vector6>::zero);
    //     Field<vector6> bPrime(x.size(), pTraits<vector6>::zero);


        
        
    //     do
    //     {
    //         rhoOld = rho;

    //         // Update search directions
    //         rho = gSumProd(rw, r);

    //         beta = rho/rhoOld*(alpha/omega);

    //         // Restart if breakdown occurs
    //         if (rho == 0)
    //         {
    //             rw = r;
    //             rho = gSumProd(rw, r);

    //             alpha = 0;
    //             omega = 0;
    //             beta = 0;
    //         }

    //         forAll (p, i)
    //         {
    //             p[i] = r[i] + beta*p[i] - beta*omega*v[i];
    //         }

    //         preconPtr_->precondition(ph, p);

    //         matrix.Amul(v, ph);
            
    //         alpha = rho/gSumProd(rw, v);

    //         forAll (s, i)
    //         {
    //             s[i] = r[i] - alpha*v[i];
    //         }

    //         // Bug fix, Alexander Monakov, 11/Jul/2012
    //         preconPtr_->precondition(sh, s);
    //         matrix.Amul(t, sh);
    //         omega = gSumProd(t, s)/gSumProd(t, t);

    //         forAll (x, i)
    //         {
    //             x[i] = x[i] + alpha*ph[i] + omega*sh[i];
    //         }

    //         forAll (r, i)
    //         {
    //             r[i] = s[i] - omega*t[i];
    //         }

    //         solverPerf.finalResidual() = cmptDivide(gSum(cmptMag(r)), norm);
    //         // solverPerf.finalResidual() = gSum(cmptMag(r))/norm;
    //         solverPerf.nIterations()++;
    //     } while (!this->stop(solverPerf));
    // }

    
    // Calculate final residual
    {
        Field<vector6> p(x.size());
        matrix.Amul(p, x);
        Field<vector6> r(b - p);

        solverPerf.finalResidual() = cmptDivide(gSum(cmptMag(r)),norm);
        solverPerf.nIterations()++;
    }
    
    return solverPerf;
}


// ************************************************************************* //
