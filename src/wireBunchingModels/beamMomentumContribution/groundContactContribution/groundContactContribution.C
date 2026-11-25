/*---------------------------------------------------------------------------*\
License
    This file is part of beamFoam.

    beamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    beamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with beamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "groundContactContribution.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(groundContactContribution, 0);
addToRunTimeSelectionTable
(
    beamMomentumContribution, groundContactContribution, dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


groundContactContribution::groundContactContribution
(
    const word& name,
    const dictionary& dict
)
:
    beamMomentumContribution
    (
        name,
        dict
    ),
    beamMomentumContribDict_(dict.subDict(name + "Coeffs")),
    kNormal_
    (
        readScalar(beamMomentumContribDict_.lookup("kNormal"))
    ),
    cNormal_
    (
        readScalar(beamMomentumContribDict_.lookup("cNormal"))
    ),
    kTangential_
    (
        readScalar(beamMomentumContribDict_.lookup("kTangential"))
    ),
    muFriction_
    (
        readScalar(beamMomentumContribDict_.lookup("muFriction"))
    ),
    groundPlane_(readGroundPlane(beamMomentumContribDict_))
{
    Info<< "Found beamMomentumContribution type: " << typeName << endl;
}


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

plane groundContactContribution::readGroundPlane(const dictionary& dict)
{
    const bool pointAndNormal =
        dict.found("groundPoint") || dict.found("groundNormal");

    if (dict.found("groundZ") == pointAndNormal)
    {
        FatalIOErrorInFunction(dict)
            << "Specify the ground plane either with groundZ, or with"
            << " groundPoint and groundNormal" << exit(FatalIOError);
    }

    if (pointAndNormal)
    {
        // The plane normalises the normal, and stops if it is zero
        return plane
        (
            dict.get<point>("groundPoint"),
            dict.get<vector>("groundNormal")
        );
    }

    // Short form: the horizontal plane z = groundZ
    return plane
    (
        point(0, 0, readScalar(dict.lookup("groundZ"))),
        vector(0, 0, 1)
    );
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<vectorField> groundContactContribution::cellPositions
(
    const beamModel& bm
) const
{
    const volVectorField& W = bm.solutionW();
    const fvMesh& mesh = W.mesh();

    // The mesh cell-centres plus the reference displacement (registered by
    // the beam model) and the displacement
    return tmp<vectorField>::New
    (
        mesh.C().internalField()
      + mesh.lookupObject<volVectorField>("refW").internalField()
      + W.internalField()
    );
}


scalar groundContactContribution::normalForce
(
    const scalar R,
    const scalar d,
    const scalar Un
) const
{
    // The spring pushes the penetrating cell out, and the damper resists its
    // rebound (velocity along the normal). The seabed can only push, so a
    // damper stronger than the spring does not pull the cell in
    return max
    (
        (2*R*kNormal_*(-d))
      - (2*R*cNormal_*max(Un, 0)),
        0
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<Field<scalarSquareMatrix>> groundContactContribution::diagCoeff
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
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
    Field<scalarSquareMatrix>& result = tresult.ref();

    // The diagonal holds the derivative of the residual (internal plus
    // external force minus inertia) with respect to W, like the implicit
    // inertia coefficient. On a cell in contact, the normal contact force is
    // normalForce*L along the ground normal n, and its spring part
    // 2*R*kNormal*(-d)*L depends on W through the signed distance d, with
    // dd/dW = n. Its derivative is therefore -2*R*kNormal*L*(n n). It is
    // added only where the clipped normal force is positive: where the
    // damper has cut the force to zero, the force does not depend on W. The
    // damping and friction terms stay explicit
    const vectorField coord(cellPositions(bm));
    const vector& n = groundPlane_.normal();
    const scalar R = bm.R();
    const volScalarField& L = bm.L();

    forAll(result, cellI)
    {
        const scalar d = groundPlane_.signedDistance(coord[cellI]);

        if (d < 0 && normalForce(R, d, U[cellI] & n) > 0)
        {
            const scalar k = 2*R*kNormal_*L[cellI];

            for (direction i = 0; i < vector::nComponents; ++i)
            {
                for (direction j = 0; j < vector::nComponents; ++j)
                {
                    result[cellI](i, j) -= k*n[i]*n[j];
                }
            }
        }
    }

    return tresult;
}


tmp<vectorField> groundContactContribution::linearMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
)
{
    const volVectorField& W = bm.solutionW();
    // Take a reference to the mesh
    const fvMesh& mesh = W.mesh();

    // Current position of the beam CV cell-centres
    const vectorField coord(cellPositions(bm));

    // Unit normal of the ground plane, pointing from the ground into the fluid
    const vector& n = groundPlane_.normal();

    // Beam Radius
    const scalar R = bm.R();

    // Prepare the result
    tmp<vectorField> tresult(new vectorField(mesh.nCells(), vector::zero));
    vectorField& result = tresult.ref();

    // TO-DO: The friction below only resists sliding along the beam, not
    // sliding across it on the ground

    // Unit tangents to beam centreline at beam CV cell-centres
    const vectorField dRdScellHat(unitTangents(bm));

    // Tangential component of velocity vector
    vectorField Ut
    (
        (
            (U.internalField() & dRdScellHat)
            *dRdScellHat
        )
    );

    vectorField UtHat (Ut/(mag(Ut) + SMALL));

    // Beam cell lengths: the forces below are per unit length
    const volScalarField& L = bm.L();

    // Initialise beam cells in contact with ground
    label cellsInContact = 0;

    forAll(result, cellI)
    {
        // Signed distance of the cell-centre from the ground plane,
        // negative below the ground
        const scalar d = groundPlane_.signedDistance(coord[cellI]);

        if (d < 0)
        {
            cellsInContact++;
            // Seabed force per unit length, along the ground normal
            const scalar f_gc_normal = normalForce(R, d, U[cellI] & n);

            vector f_gc_tangential(vector::zero);

            if
            (
                kTangential_*2.0*R*mag(Ut[cellI])
             >= muFriction_*f_gc_normal
            )
            {
                // Max value of friction as per Coulomb's law
                f_gc_tangential =
                    -muFriction_*f_gc_normal*UtHat[cellI];
            }
            else
            {
                f_gc_tangential = -2.0*R*kTangential_*Ut[cellI];
            }

            // The solver adds external forces to the source with a minus
            // sign, as for gravity, so subtract the seabed force. The forces
            // above are per unit length, so multiply by the cell length, as
            // for the distributed load q and the Morison drag
            result[cellI] -= (f_gc_normal*n + f_gc_tangential)*L[cellI];
        }
     }

    Info<< "Number of cells in contact : " << cellsInContact << endl;

    return tresult;
}


tmp<vectorField> groundContactContribution::angularMomentumSource
(
    const beamModel& bm,
    const volVectorField& U,
    const volVectorField& Accl
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
