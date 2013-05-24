/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::TriSurfaceMeshDistanceCalculator

Description

SourceFiles
    TriSurfaceMeshDistanceCalculator.C
  
Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modeling And Analysis Group
    Center of Smart Interfaces
    TU Darmstadt

\*---------------------------------------------------------------------------*/

#ifndef TriSurfaceMeshDistanceCalculator_H
#define TriSurfaceMeshDistanceCalculator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class TriSurfaceMeshDistanceCalculator Declaration
\*---------------------------------------------------------------------------*/

class TriSurfaceMeshDistanceCalculator
{
    //- List of indices to the nearest elements to cells.
    DynamicList<pointIndexHit> cellsElementNearest_; 

    //- List of indices to the nearest elements to points.
    DynamicList<pointIndexHit> pointsElementNearest_; 

    //- Cell search distance field used by the octree search algorithm.
    autoPtr<volScalarField> cellSearchDistPtr_;

    //- Point search distance field used by the octree search algorithm.
    autoPtr<pointScalarField> pointSearchDistPtr_;

    //- Initialize the cell search distance field. 
    void initCellSearchDistance(const fvMesh& mesh);

    //- Initialize the point search distance field. 
    void initPointSearchDistance(const fvMesh& mesh);

    // Narrow band width in cell layers.
    label bandwidth_; 

    // Compute the width of the cell search stencil for the octree search. 
    void calcCellSearchDistance(const fvMesh& mesh); 
    // Compute the width of the point search stencil for the octree search. 
    void calcPointSearchDistance(const fvMesh& mesh); 

public:

    TriSurfaceMeshDistanceCalculator(label bandwidth);

    // Member Functions
    
    template <typename Connection, typename NarrowBandPropagation>
    void calcCentresToElementsDistance
    (
        volScalarField& Psi, 
        const Connection& connection, 
        NarrowBandPropagation enforceNarrowBand 
    ); 

    template <typename Connection>
    void calcPointsToElementsDistance
    (
        scalarField& psi, 
        const Connection& connection
    ); 

    void setBandwidth(label bandwidth)
    {
        bandwidth_ = bandwidth; 
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 
} // End namespace Foam

#ifdef NoRepository
#   include "TriSurfaceMeshDistanceCalculator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
