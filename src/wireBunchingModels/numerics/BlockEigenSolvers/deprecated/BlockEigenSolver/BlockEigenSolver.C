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

// #include "blockLduSolvers.H"
#include "BlockEigenSolver.H"
#include "OFstream.H"
#include "polyMesh.H"
#include <unsupported/Eigen/SparseExtra>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BlockEigenSolver, 0);
}

// * * * * * * * * * * * Protected Data Functions * * * * * * * * * * * * * //

bool Foam::BlockEigenSolver::checkTwoD() const
{
    // // Check if mesh is 2-D or 3-D
    // // This should be OK for multi-region cases as all regions should be 2-D in
    // // same direction
    // const polyMesh* meshPtr = NULL;
    // if (matrix_.mesh().thisDb().parent().foundObject<polyMesh>("solid"))
    // {
    //     meshPtr =
    //         &matrix_.mesh().thisDb().parent().lookupObject<polyMesh>("solid");
    // }
    // else
    // {
    //     meshPtr =
    //         &matrix_.mesh().thisDb().parent().lookupObject<polyMesh>("region0");
    // }
    // const polyMesh& mesh = *meshPtr;

    // const Vector<label>& geomD = mesh.geometricD();

    // label nDir = 0;
    // forAll(geomD, dirI)
    // {
    //     if (geomD[dirI] > SMALL)
    //     {
    //         nDir++;
    //     }
    // }

    // bool twoD = false;
    // if (nDir == 2)
    // {
    //     if (mesh.solutionD()[vector::Z] > -1)
    //     {
    //         FatalErrorIn
    //         (
    //             "Foam::BlockSolverPerformance<Foam::vector>"
    //             "Foam::BlockEigenSolver::solve"
    //             "("
    //             "    Field<Foam::vector>& U,"
    //             "    const Field<Foam::vector>& blockB"
    //             ")"
    //         )   << this->typeName << ": for 2-D models, the empty direction "
    //             << "must be z!" << abort(FatalError);
    //     }

    //     twoD = true;
    // }
    // else if (nDir != 3)
    // {
    //     FatalErrorIn(this->typeName + " solve")
    //         << "solver only implemented for 2-D and 3-D models"
    //         << abort(FatalError);
    // }

    return false;
}


int Foam::BlockEigenSolver::calcDegreesOfFreedom
(
    const BlockLduMatrix<vector6>& matrix,
    const bool twoD
) const
{
    // Number of vectors on the diagonal
    const int diagSize = matrix.diag().asSquare().size();

    if (twoD)
    {
        // Two components in each vector
        return 4*diagSize;
    }

    // Info << "diagSize: " << diagSize << endl;
    // Three components in each vector
    return 6*diagSize;
}


void Foam::BlockEigenSolver::convertFoamMatrixToEigenMatrix
(
    const BlockLduMatrix<vector6>& matrix,
    Eigen::SparseMatrix<scalar>& A
)
{
    if (BlockLduSolver::debug)
    {
        Info<< this->typeName
            << ": copying matrix coefficients into Eigen format"
            << endl;
    }

    // Matrix addressing
    const unallocLabelList& lowerAddr = matrix.mesh().lduAddr().lowerAddr();
    const unallocLabelList& upperAddr = matrix.mesh().lduAddr().upperAddr();

    // Grab matrix diagonal and off-diagonals
    const Field<tensor6>& d = matrix.diag().asSquare();
    const Field<tensor6>& l = matrix.lower().asSquare();
    const Field<tensor6>& u = matrix.upper().asSquare();

    // Check if the matrix is from a 2-D case
    const bool twoD = checkTwoD();
    // Info << twoD << endl;

    // Calculate te number of degrees of freedom
    const int m = calcDegreesOfFreedom(matrix, twoD);

    // Info << "Number of degrees of freedom: " << m << endl;
    
    // Create coefficient matrix: we must copy coeffs from ldu storage
    // This can be costly
    std::vector< Eigen::Triplet<scalar> > coefficients;
    coefficients.reserve(m);

    // Insert coefficients

    // Insert diagonal
    label ID = -1;
    forAll(d, dI)
    {
        if (twoD)
        {
            // 2-D
            ID = 4*dI;
        }
        else
        {
            // 3-D
            ID = 6*dI;
        }

        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID, d[dI](0,0)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID+1, d[dI](0,1)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID+2, d[dI](0,2)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID+3, d[dI](0,3)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID+4, d[dI](0,4)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID, ID+5, d[dI](0,5)));

        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID, d[dI](1,0)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID+1, d[dI](1,1)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID+2, d[dI](1,2)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID+3, d[dI](1,3)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID+4, d[dI](1,4)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+1, ID+5, d[dI](1,5)));

        coefficients.push_back(Eigen::Triplet<scalar>(ID+2, ID, d[dI](2,0)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+2, ID+1, d[dI](2,1)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+2, ID+2, d[dI](2,2)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+2, ID+3, d[dI](2,3)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+2, ID+4, d[dI](2,4)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+2, ID+5, d[dI](2,5)));

        coefficients.push_back(Eigen::Triplet<scalar>(ID+3, ID, d[dI](3,0)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+3, ID+1, d[dI](3,1)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+3, ID+2, d[dI](3,2)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+3, ID+3, d[dI](3,3)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+3, ID+4, d[dI](3,4)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+3, ID+5, d[dI](3,5)));

        coefficients.push_back(Eigen::Triplet<scalar>(ID+4, ID, d[dI](4,0)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+4, ID+1, d[dI](4,1)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+4, ID+2, d[dI](4,2)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+4, ID+3, d[dI](4,3)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+4, ID+4, d[dI](4,4)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+4, ID+5, d[dI](4,5)));

        coefficients.push_back(Eigen::Triplet<scalar>(ID+5, ID, d[dI](5,0)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+5, ID+1, d[dI](5,1)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+5, ID+2, d[dI](5,2)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+5, ID+3, d[dI](5,3)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+5, ID+4, d[dI](5,4)));
        coefficients.push_back(Eigen::Triplet<scalar>(ID+5, ID+5, d[dI](5,5)));        
    }

    // Insert off-diagonal
    label rowI = -1;
    label columnI = -1;
    forAll(u, uI)
    {
        const label own = lowerAddr[uI];
        const label nei = upperAddr[uI];

        const tensor6& upper = u[uI];
        const tensor6& lower = l[uI];

        if (twoD)
        {
            rowI = 4*own;
            columnI = 4*nei;
        }
        else
        {
            rowI = 6*own;
            columnI = 6*nei;
        }

        // Upper
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI, upper(0,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI+1, upper(0,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI+2, upper(0,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI+3, upper(0,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI+4, upper(0,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI, columnI+5, upper(0,5))
        );

        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI, upper(1,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI+1, upper(1,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI+2, upper(1,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI+3, upper(1,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI+4, upper(1,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+1, columnI+5, upper(1,5))
        );

        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+2, columnI, upper(2,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+2, columnI+1, upper(2,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+2, columnI+2, upper(2,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+2, columnI+3, upper(2,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+2, columnI+4, upper(2,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+2, columnI+5, upper(2,5))
        );
        
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+3, columnI, upper(3,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+3, columnI+1, upper(3,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+3, columnI+2, upper(3,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+3, columnI+3, upper(3,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+3, columnI+4, upper(3,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+3, columnI+5, upper(3,5))
        );
        
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+4, columnI, upper(4,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+4, columnI+1, upper(4,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+4, columnI+2, upper(4,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+4, columnI+3, upper(4,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+4, columnI+4, upper(4,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+4, columnI+5, upper(4,5))
        );

        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+5, columnI, upper(5,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+5, columnI+1, upper(5,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+5, columnI+2, upper(5,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+5, columnI+3, upper(5,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+5, columnI+4, upper(5,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(rowI+5, columnI+5, upper(5,5))
        );


        // Lower
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI, lower(0,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI+1, lower(0,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI+2, lower(0,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI+3, lower(0,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI+4, lower(0,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI, rowI+5, lower(0,5))
        );

        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI, lower(1,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI+1, lower(1,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI+2, lower(1,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI+3, lower(1,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI+4, lower(1,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+1, rowI+5, lower(1,5))
        );

        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+2, rowI, lower(2,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+2, rowI+1, lower(2,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+2, rowI+2, lower(2,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+2, rowI+3, lower(2,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+2, rowI+4, lower(2,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+2, rowI+5, lower(2,5))
        );
        
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+3, rowI, lower(3,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+3, rowI+1, lower(3,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+3, rowI+2, lower(3,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+3, rowI+3, lower(3,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+3, rowI+4, lower(3,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+3, rowI+5, lower(3,5))
        );
        
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+4, rowI, lower(4,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+4, rowI+1, lower(4,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+4, rowI+2, lower(4,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+4, rowI+3, lower(4,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+4, rowI+4, lower(4,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+4, rowI+5, lower(4,5))
        );

        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+5, rowI, lower(5,0))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+5, rowI+1, lower(5,1))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+5, rowI+2, lower(5,2))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+5, rowI+3, lower(5,3))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+5, rowI+4, lower(5,4))
        );
        coefficients.push_back
        (
            Eigen::Triplet<scalar>(columnI+5, rowI+5, lower(5,5))
        );
    }

    // Insert triplets into the matrix
    A.setFromTriplets(coefficients.begin(), coefficients.end());
}

void Foam::BlockEigenSolver::writeLinearSystemToMatlabFiles
(
    const Eigen::SparseMatrix<scalar>& A,
    const Eigen::Matrix<scalar, Eigen::Dynamic, 1>& b
) const
{
    if (writeMatlabFiles_)
    {
        Info<< type() << ": writing sparse matrix to "
            << "matlabSparseMatrix.txt and source to matlabSource.txt"
            << nl
            << "The system can be read and solved in Matlab or Octave with "
            << "commands:" << nl
            << nl
            << "    load matlabSparseMatrix.txt;" << nl
            << "    A = spconvert(matlabSparseMatrix);" << nl
            << "    B = dlmread('matlabSource.txt', ' ');" << nl
            << "    x = A\\B;" << nl
            << endl;

        // Write matrix
        Eigen::saveMarket(A, "matlabSparseMatrix.txt");

        // Write source
        OFstream sourceFile("matlabSource.txt");
        for (int rowI = 0; rowI < A.rows(); rowI++)
        {
            sourceFile
                << b(rowI) << endl;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::BlockEigenSolver::BlockEigenSolver
(
    const word& fieldName,
    const BlockLduMatrix<vector6>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<vector6>(fieldName, matrix, dict),
    writeMatlabFiles_
    (
        dict.lookupOrDefault<Switch>("writeMatlabFiles", false)
    )
{
    // Info<< type() << " : writeMatlabFiles is " << writeMatlabFiles_ << endl;
}


// ************************************************************************* //
