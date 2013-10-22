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
    Foam::TriSurfaceMeshCalculator

Description
    Calculates front<->mesh properties based on the octree-connectivity 
    between the front and the mesh provided by the triSurfaceMesh class.

SourceFiles
    TriSurfaceMeshCalculator.C
  
Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modeling And Analysis Group
    Center of Smart Interfaces
    TU Darmstadt

\*---------------------------------------------------------------------------*/

#ifndef TriSurfaceMeshCalculator_H
#define TriSurfaceMeshCalculator_H

#include "DynamicField.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class TriSurfaceMeshCalculator Declaration
\*---------------------------------------------------------------------------*/

class TriSurfaceMeshCalculator
{
    // FIXME: add autPtr for the frontMesh. 
    //- Front Mesh autoPtr
    mutable autoPtr<triSurfaceMesh> frontMeshPtr_; 

    //- List of indices to the nearest elements to cells.
    DynamicList<pointIndexHit> cellsElementNearest_; 

    //- List of indices to the nearest elements to points.
    DynamicList<pointIndexHit> pointsElementNearest_; 

    //- Cell search distance field used by the octree search algorithm.
    autoPtr<volScalarField> cellSearchDistSqrPtr_;

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

    // Loose point in cell query: if the point is on a face, the function returns true. 
    bool pointInCell(const point& p, label cellI, const fvMesh& mesh) const; 

    // Heuristic cell locator, based on a seed label.  
    label findCell (const point&p, label cellI, const fvMesh& mesh) const; 

    // TODO: recode these, force deep copy and not ownership transfer using new autoPtr 
    TriSurfaceMeshCalculator(const TriSurfaceMeshCalculator&); 
    void operator=(const TriSurfaceMeshCalculator&); 


public:

    TriSurfaceMeshCalculator(label bandwidth);

    TriSurfaceMeshCalculator(const fvMesh& mesh, label bandWidth); 

    // Member Functions
    
    template <typename NarrowBandPropagation>
    void calcCentresToElementsDistance
    (
        volScalarField& Psi, 
        const triSurfaceFront& front, 
        NarrowBandPropagation enforceNarrowBand
    ); 

    template<typename Mesh, typename NarrowBandPropagation>
    void calcPointsToElementsDistance
    (
        pointScalarField& psi, 
        const triSurfaceFront& front, 
        const Mesh& mesh, 
        NarrowBandPropagation enforceNarrowBand
    ); 

    void calcFrontVelocity
    (
        DynamicField<vector>& frontVelocity, 
        triSurfaceFront& front, 
        const volVectorField& U
    );

    void setBandwidth(label bandwidth)
    {
        bandwidth_ = bandwidth; 
    }

    const volScalarField& getCellSearchDistSqr()
    {
        return cellSearchDistSqrPtr_(); 
    }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 
} // End namespace Foam

#ifdef NoRepository
#   include "TriSurfaceMeshCalculator.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
