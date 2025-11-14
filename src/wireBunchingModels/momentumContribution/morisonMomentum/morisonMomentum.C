/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "morisonMomentum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(morisonMomentum, 0);
addToRunTimeSelectionTable(momentumContribution, morisonMomentum, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


morisonMomentum::morisonMomentum
(
    const Time& runTime
)
:
    momentumContribution(runTime) //,
    // Um_
    // (
    //     this->subDict(typeName + "Coeffs").lookup("Um")
    // ),
    // Tsoft_
    // (
    //     readScalar(this->subDict(typeName + "Coeffs").lookup("Tsoft"))
    // )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<Field<scalarSquareMatrix>> morisonMomentum::diagCoeff
(
    const volVectorField& W, const volVectorField& Theta
)
{
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<Field<scalarSquareMatrix>> tresult
     (
         new Field<scalarSquareMatrix>
         (
             mesh.nCells(), scalarSquareMatrix(6, 6, 0.0)
         )
     );
    // Field<scalarSquareMatrix>& result = tresult.ref();

    return tresult;
}


tmp<vectorField> morisonMomentum::linearMomentumSource
(
    const volVectorField& W, const volVectorField& Theta
)
{
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    // vectorField& result = tresult.ref();

    return tresult;
}


tmp<vectorField> morisonMomentum::angularMomentumSource
(
    const volVectorField& W, const volVectorField& Theta
)
{
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    // vectorField& result = tresult.ref();

    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
