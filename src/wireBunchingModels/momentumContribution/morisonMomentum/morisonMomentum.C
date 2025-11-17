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
    momentumContribution(runTime),
    Cdn_
    (
        readScalar(this->subDict(typeName + "Coeffs").lookup("Cdn"))
    ),
    Cdt_
    (
        readScalar(this->subDict(typeName + "Coeffs").lookup("Cdt"))
    )
{
    Info<< "Found momentumContribution type: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<Field<scalarSquareMatrix>> morisonMomentum::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U
    // const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Prepare the result
    tmp<Field<scalarSquareMatrix>> tresult
     (
        new Field<scalarSquareMatrix>
        (
            mesh.nCells(), scalarSquareMatrix(6, 0.0)
        )
     );
    // Field<scalarSquareMatrix>& result = tresult.ref();

    return tresult;
}


tmp<vectorField> morisonMomentum::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U
    // const volVectorField& Accl
)
{
    // Take a reference to the mesh
    const fvMesh& mesh = bm.solutionW().mesh();

    // Flag to check the ddtSchemeName
    const word ddtSchemeName
    (
        mesh.ddtSchemes().found("ddt(W)")
        ?
        (mesh.ddtSchemes().lookup("ddt(W)"))
        :
        (mesh.ddtSchemes().lookup("default"))
    );

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    vectorField& result = tresult.ref();

    //- Drag forces due to Morison's Equation
    if (ddtSchemeName != "steadyState")
    {
        // Create spline using current beam points and tangents data
        HermiteSpline spline
        (
            bm.currentBeamPoints(),
            bm.currentBeamTangents()
        );

        // Evaluate dRdS - tangents to beam centreline at beam CV cell-centres
        const vectorField& dRdScell = spline.midPointDerivatives();

        // Tangential component of velocity vector
        vectorField Ut
        (
            (
                (U.internalField() & dRdScell)
                *dRdScell
            )
        );

        vectorField UtHat (Ut/(mag(Ut) + SMALL));

        // Normal component of velocity vector
        vectorField Un
        (
            (
                U.internalField()
                - (
                    (U.internalField() & dRdScell)
                    *dRdScell
                )
            )
        );

        vectorField UnHat (Un/(mag(Un) + SMALL));

        const scalar rho = bm.rho().value();
        const scalar R = bm.R();
        const volScalarField& L = bm.L();

        // Scalar values of drag force (normal and tangential)
        const scalarField Fdn(rho*Cdn_*R*L*(Un & Un));
        const scalarField Fdt(rho*Cdt_*R*L*(Ut & Ut));

    //     // Explicit drag forces included in the source vector
    //     forAll(source, cellI)
    //     {
    //         source[cellI](0,0) += Fdn[cellI]*UnHat[cellI].component(0);
    //         source[cellI](1,0) += Fdn[cellI]*UnHat[cellI].component(1);
    //         source[cellI](2,0) += Fdn[cellI]*UnHat[cellI].component(2);

    //         source[cellI](0,0) += Fdt[cellI]*UtHat[cellI].component(0);
    //         source[cellI](1,0) += Fdt[cellI]*UtHat[cellI].component(1);
    //         source[cellI](2,0) += Fdt[cellI]*UtHat[cellI].component(2);
    //     }

        forAll(result, cellI)
        {
            result[cellI] =
                Fdn[cellI]*UnHat[cellI]
              + Fdt[cellI]*UtHat[cellI];
        }
    }
    else
    {
        FatalErrorInFunction
        << "Momentum contribution due to Morison drag - type: "
        << typeName
        << " for time scheme " << ddtSchemeName
        << " is not defined!!\n"
        << "Set d2dt2Scheme and ddtScheme in system/fvSchemes as "
        << "'Euler' or 'Newmark'"
        << "to include drag force contributions"
        << abort(FatalError);
    }

    return tresult;
}


tmp<vectorField> morisonMomentum::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U
    // const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
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
