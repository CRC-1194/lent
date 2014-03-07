/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::frontVelocityCalculator

Description
    Interface for the front velocity calculation. 

SourceFiles
    frontVelocityCalculator.C

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#ifndef frontVelocityCalculator_H
#define frontVelocityCalculator_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"
#include "triSurfaceFrontFields.H"
#include "lentMeshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontVelocityCalculator Declaration
\*---------------------------------------------------------------------------*/

class frontVelocityCalculator
    :
        public refCount
{

    //tmp<lentMeshSearch> searchAlgTmp_; 

public:

    TypeName ("frontVelocityCalculator"); 

    declareRunTimeSelectionTable (
        tmp,  
        frontVelocityCalculator,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors

        frontVelocityCalculator(const dictionary& configDict);

    // Selectors

        static tmp<frontVelocityCalculator> New(const dictionary& configDict);

    // Destructor
        virtual ~frontVelocityCalculator();

    // Member Functions

        //bool pointIsInCell(
            //const point p, 
            //const label cellLabel, 
            //const fvMesh& mesh, 
            //scalar tolerance=SMALL
        //) const 
        //{
            //return searchAlgTmp_->pointIsInCell(p, cellLabel, mesh, tolerance); 
        //} 

        //label cellContainingPoint(
            //const point& p, 
            //const fvMesh& mesh, 
            //const label seedCell = -1
        //) const
        //{
            //return searchAlgTmp_->cellContainingPoint(p, mesh, seedCell);  
        //} 

        //virtual void updateMeshCells(
            //labelList& meshCells, 
            //const triSurfaceFront& front,
            //const fvMesh& mesh
        //);
        
        virtual void calcFrontVelocity(
            triSurfaceFrontVectorField& frontVelocity, 
            const volVectorField& U,
            labelList& elementCells 
        ) const = 0;  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
