/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    circle

Description
    Circle cross-section model

SourceFiles
    circle.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef circle_H
#define circle_H

#include "crossSectionModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class circle Declaration
\*---------------------------------------------------------------------------*/

class circle
:
    public crossSectionModel
{
    // Private data

        //- Cross-section model dictionary
        dictionary crossSectionModelDict_;
    
        //- Circle radius
        scalar radius_;
	
	//- Area scaling factor
	scalar scalingArea_;
	
	//- Stiffness (or Moment of Inertia) scaling factor
	// Since E cannot be changed for beams this scaling
	// factor is artificially added
	scalar scalingMI_;

        //- Local mesh points
        vectorField points_;
    
    // Private Member Functions

        //- Disallow copy construct
        circle(const circle&);

        //- Disallow default bitwise assignment
        void operator=(const circle&);

public:

    //- Runtime type information
    TypeName("circle");


    // Constructors

        //- Construct from dictionary
        circle
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~circle()
        {}


    // Member Functions

        //- Return cross-section area
        virtual scalar radius(const scalar theta) const
        {
            return radius_;
        }
    
        //- Return cross-section area
        virtual scalar R() const
        {
            return radius_;
        }
    
        //- Return cross-section area
        virtual scalar A() const
        {
            return scalingArea_*M_PI*sqr(radius_);
        }
    
        //- Return cross-section moment of inertia
        virtual scalar Ixx() const
        {
            return scalingMI_*M_PI*pow(radius_, 4)/4;
        }
    
        //- Return cross-section moment of inertia
        virtual scalar Iyy() const
        {
            return scalingMI_*M_PI*pow(radius_, 4)/4;
        }
    
        //- Return cross-section Saint-Venant torsion constant
        virtual scalar IT() const
        {
            return scalingMI_*M_PI*pow(radius_, 4)/2;
        }
    
        //- Return cross-section coefficient of proportionality
        //  between fully plastic axial force and fully plastic moment
        virtual scalar CNM() const
        {
            return 3*M_PI/(4*radius_);
        }

        //- Circumferential points
        virtual tmp<vectorField> circumferentialPoints() const
        {
            label nTheta = 32;
            
            tmp<vectorField> tCircumPoints
            (
                new vectorField(nTheta, vector::zero)
            );

            scalar dTheta = -2*M_PI/nTheta;

            forAll(tCircumPoints(), pointI)
            {
                tCircumPoints()[pointI].x() = radius_*::cos(pointI*dTheta);
                tCircumPoints()[pointI].y() = radius_*::sin(pointI*dTheta);
                tCircumPoints()[pointI].z() = 0;
            }

            return tCircumPoints;
        }

        //- Number of circumferential points
        virtual label nCircumferentialPoints () const
        {
            label nTheta = 32;
            return nTheta;
        }

        //- Structured mesh points
        virtual const vectorField& points() const
        {
            return points_;
        }

        //- Return Green-Lagrangian strain at structured mesh points
        //  for provided average tranlational and rotational strains
        virtual tmp<vectorField> greenLagrangianStrain
        (
            const vector& Gamma,
            const vector& K
        ) const;
    
        //- Return resultant cross-section force
        virtual tmp<scalarField> resultantTotalForce
        (
            const vectorField& S
        ) const
        {
            tmp<scalarField> trF(new scalarField(6, 0));
            return trF;
        }
    
        //- Return resultant cross-section force
        virtual vector resultantForce(const vectorField& S) const
        {
            vector rF = vector::zero;
            return rF;
        }
    
        //- Return resultant cross-section moment
        virtual vector resultantMoment(const vectorField& S) const
        {
            vector rM = vector::zero;
            return rM;
        }            

        //- Return cross-section resultant tangent matrices
        void resultantTangentMatrices
        (
            // const vectorField& S,
            const diagTensorField& dSdE,
            tensor& DQDGamma,
            tensor& DQDK,
            tensor& DMDGamma,
            tensor& DMDK
        ) const {}
    
        //- Write cross-section scalar field in VTK format
        virtual void writeVTK
        (
            const fileName& fn,
            const scalarField& sf0,
            const scalarField& sf1,
            const vectorField& vf
        ) const {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
