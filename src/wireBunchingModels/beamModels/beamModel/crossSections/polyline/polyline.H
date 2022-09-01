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
    polyline

Description
    Polyline cross-section model

SourceFiles
    polyline.C

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef polyline_H
#define polyline_H

#include "crossSectionModel.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "faceList.H"
#include "quadrilateral.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class polyline Declaration
\*---------------------------------------------------------------------------*/

class polyline
:
    public crossSectionModel
{
    // Private data

        //- Cross-section model dictionary
        dictionary crossSectionModelDict_;
    
        //- Circumferential points
        vectorField circumferentialPoints_;

        //- Minimal cross-section radius
        scalar minRadius_;
    
        //- Cross-section area
        scalar A_;

        //- Cross-section second moment of area
        scalar Ixx_;
    
        //- Cross-section second moment of area
        scalar Iyy_;
    
        //- Number of circumferential points
        label nCircumferentialPoints_;
    
        //- Number of segments per side
        label nSegmentsPerSide_;

        //- Number of segments in radial direction
        label nRadialSegments_;
    
        //- Local mesh points
        vectorField points_;

        //- Local mesh faces
        // faceList quadFaces_;
        List<quadrilateral> quadFaces_;
        // faceList faces_;
    
    // Private Member Functions

        //- Calculate cross-section area
        scalar area() const;
    
        //- Calculate second moment of area
        void calcSecondMomentOfArea(scalar& Ixx, scalar& Iyy);
    
        //- Calculate second moment of area
        scalar IYY() const;
    
        //- Disallow copy construct
        polyline(const polyline&);

        //- Disallow default bitwise assignment
        void operator=(const polyline&);

public:

    //- Runtime type information
    TypeName("polyline");


    // Constructors

        //- Construct from dictionary
        polyline
        (
            const word& name,
            const dictionary& dict
        );

    // Destructor

        virtual ~polyline()
        {}


    // Member Functions

        //- Definitian of the coross-section
        //  boundry line in polar coordinates
        virtual scalar radius(const scalar theta) const;
    
        //- Return cross-section radius
        virtual scalar R() const
        {
            return minRadius_;
        }
    
        //- Return cross-section area
        virtual scalar A() const
        {
            // return M_PI*sqr(minRadius_);
            return A_;
        }
    
        //- Return cross-section moment of inertia
        virtual scalar Ixx() const
        {
            // return M_PI*pow(minRadius_, 4)/4;
            return Ixx_;
        }

        //- Return cross-section moment of inertia
        virtual scalar Iyy() const
        {
            // return M_PI*pow(minRadius_, 4)/4;
            return Iyy_;
        }
    
        //- Return cross-section Saint-Venant torsion constant
        virtual scalar IT() const
        {
            // return M_PI*pow(minRadius_, 4)/2;
            return Ixx_ + Iyy_;
        }

        //- Return cross-section coefficient of proportionality
        //  between fully plastic axial force and fully plastic moment
        virtual scalar CNM() const
        {
            return 4.0/minRadius_;
        }

        //- Circumferential points
        virtual tmp<vectorField> circumferentialPoints() const;
            
        //- Number of circumferential points
        virtual label nCircumferentialPoints () const
        {
            return nCircumferentialPoints_;
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

        //- Return cross-section resultant force
        virtual tmp<scalarField> resultantTotalForce
        (
            const vectorField& S
        ) const;
    
        //- Return cross-section resultant force
        virtual vector resultantForce(const vectorField& S) const;

        //- Return cross-section resultant moment
        virtual vector resultantMoment(const vectorField& S) const;

        //- Return cross-section resultant tangent matrices
        virtual void resultantTangentMatrices
        (
            // const vectorField& S,
            const diagTensorField& dSdE,
            tensor& DQDGamma,
            tensor& DQDK,
            tensor& DMDGamma,
            tensor& DMDK
        ) const;
    
        //- Write cross-section scalar field in VTK format
        virtual void writeVTK
        (
            const fileName& fn,
            const scalarField& sf0,
            const scalarField& sf1,
            const vectorField& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
