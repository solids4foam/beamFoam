/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "beamEnergyData.H"
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOmanip.H"
#include "beamModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(beamEnergyData, 0);
    addToRunTimeSelectionTable(functionObject, beamEnergyData, dictionary);
}
}

// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

bool Foam::functionObjects::beamEnergyData::writeData()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const beamModel& beam =
        mesh.objectRegistry::parent().lookupObject<beamModel>("beamProperties");

    // The internal strains and curvatures for calculation of internal energy
    const surfaceVectorField& Gamma =
        mesh.lookupObject<surfaceVectorField>("Gamma");
    const surfaceVectorField& K = mesh.lookupObject<surfaceVectorField>("K");
    const surfaceTensorField& CQ = mesh.lookupObject<surfaceTensorField>("CQ");
    const surfaceTensorField& CM = mesh.lookupObject<surfaceTensorField>("CM");

    // Properties of the beam like length, area, second moment of area, density
    const volScalarField& L = beam.L();
    const dimensionedScalar& rho = beam.rho();
    const dimensionedScalar& A = beam.A();
    const dimensionedScalar& Iyy = beam.Iyy();
    const dimensionedScalar& Izz = beam.Izz();
    const dimensionedScalar& J = beam.J();

    // Second moment of area scaling factor-if not specified in beamProperties
    // then the value is equal to 1.0
    const scalar kCI = beam.kCI();

    // Inertia tensor
    tensor CI
    (
        kCI*J.value(), 0, 0,
        0, kCI*Iyy.value(), 0,
        0, 0, kCI*Izz.value()
    );

    // Beam linear and angular velocity for calculation of kinetic energy
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const volVectorField& Omega = mesh.lookupObject<volVectorField>("Omega");

    // Accessing internal fields
    const vectorField& GammaI = Gamma.internalField();
    const tensorField& CQI = CQ.internalField();
    const vectorField& KI = K.internalField();
    const tensorField& CMI = CM.internalField();
    const vectorField& UI = U.internalField();
    const vectorField& OmegaI = Omega.internalField();

    // Energy parameters initialised
    scalar intE(0);
    scalar kinE_lin(0);
    scalar kinE_ang(0);
    scalar kinE(0);
    scalar tot_E(0);

    // Calculation of internal energy at internal faces
    forAll(GammaI, faceI)
    {
        intE += 0.5*L[0]*
        (
            (GammaI[faceI] & (CQI[faceI] & GammaI[faceI]))
          + (KI[faceI] & (CMI[faceI] & KI[faceI]))
        );
    }

    // Calculation of kinetic energy at cell centres
    // using cell-centre velocity fields
    forAll(UI, cellI)
    {
        kinE_lin += 0.5*rho.value()*L[0]*
        (
            A.value()*(UI[cellI] & UI[cellI])
        );
        kinE_ang += 0.5*rho.value()*L[0]*
        (
            (OmegaI[cellI] & (CI & OmegaI[cellI]))
        );
    }

    // Adding contribution of boundary fields to the internal energy
    forAll(Gamma.boundaryField(),patchI)
    {
        const vectorField& pGamma = Gamma.boundaryField()[patchI];
        const tensorField& pCQ = CQ.boundaryField()[patchI];
        const vectorField& pK = K.boundaryField()[patchI];
        const tensorField& pCM = CM.boundaryField()[patchI];

        forAll(pGamma,faceI)
        {
            intE += 0.25*L[0]*
            (
                (pGamma[faceI] & (pCQ[faceI] & pGamma[faceI]))
                + (pK[faceI] & (pCM[faceI] & pK[faceI]))
            );
        }
    }

    // Total kinetic energy
    kinE = kinE_lin + kinE_ang;

    // Total energy
    tot_E = intE + kinE;

    if (Pstream::master())
    {
        historyFilePtr_()
        << time_.time().value() << " "
        << intE << " "
        << kinE << " "
        << kinE_lin << " "
        << kinE_ang << " "
        << tot_E
        << endl;
    }

    return true;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::beamEnergyData::beamEnergyData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    name_(name),
    time_(runTime),
    regionName_(polyMesh::defaultRegion),
    historyFilePtr_()
{
    Info<< "Creating " << this->name() << " function object!" << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    // Create history file if not already created
    if (!historyFilePtr_.valid())
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(mesh.time().startTime().value());


            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"postProcessing"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            OStringStream FileName;
            FileName() << "beamEnergyData.dat";

            historyFilePtr_.reset
            (
                new OFstream(historyDir/word(FileName.str()))
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                << "Time" << " "
                << "E_{int}" << " "
                << "E_{kin}" << " "
                << "E_{kin-lin}" << " "
                << "E_{kin-ang}" << " "
                << "E_{tot}"
                << endl;
            }
        }
    }

    // read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::beamEnergyData::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::functionObjects::beamEnergyData::execute()
{
    return writeData();
}


bool Foam::functionObjects::beamEnergyData::end()
{
    return true;
}


bool Foam::functionObjects::beamEnergyData::write()
{
    return true;
}


// ************************************************************************* //
