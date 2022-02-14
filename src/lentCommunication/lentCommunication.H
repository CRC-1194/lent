/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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

Class
    Foam::lentCommunication

SourceFiles
    lentCommunication.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
        Front / Mesh communication maps.  

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/


#ifndef lentCommunication_H
#define lentCommunication_H

#include <vector>
#include <map>

#include "dictionary.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "DynamicList.H"
#include "frontMeshSearch.H"
#include "triSurfaceFront.H"
#include "fvMesh.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class lentCommunication Declaration
\*---------------------------------------------------------------------------*/

class lentCommunication
    :
        public regIOobject 
{
    // Front constant reference. 
    const triSurfaceFront& front_; 

    // Mesh constant reference. 
    const fvMesh& mesh_; 

    // KVS search algorithm used to update the maps.
    frontMeshSearch searchAlg_; 

    // Triangle -> Cell : triangle vertex in cell.
    DynamicList<label> triangleToCell_;

    // Vertex -> Cell : vertex in cell. 
    DynamicList<label> vertexToCell_;

    // Interface cell -> contained triangles (inverse of triangleToCell_)
    std::map<label, std::vector<label>> interfaceCellToTriangles_;

    // Interface cell -> contained vertices (inverse of vertexToCell_)
    std::map<label, std::vector<label>> interfaceCellToVertices_;

    //- Cell -> Nearest Triangle.  
    DynamicList<pointIndexHit> cellsTriangleNearest_;

    //- Point -> Nearest Triangle.  
    DynamicList<pointIndexHit> pointsTriangleNearest_;

public:

    static word registeredName(const triSurfaceFront& front, const polyMesh& mesh); 

    TypeName("dynamicListCommunication");

    declareRunTimeSelectionTable (
        autoPtr,
        lentCommunication,
        FrontMesh,
        (
            const triSurfaceFront& front,  
            const fvMesh& mesh
        ),
        (front, mesh)
    )

    // Constructors
    lentCommunication(
        const triSurfaceFront& front, 
        const fvMesh& mesh
    );

    lentCommunication(const lentCommunication&) = default; 

    // Selectors
    static autoPtr<lentCommunication> New(
        const dictionary& configDict, 
        const triSurfaceFront& front, 
        const fvMesh& mesh
    );

    // Destructor
    virtual ~lentCommunication() = default;

    // Member Functions
    const DynamicList<label>& triangleToCell() const 
    {
        return triangleToCell_; 
    }

    DynamicList<label>& triangleToCell() // TODO: Remove after removing mesh search from the velocity interpolator. Follow through in lentMethod class. TM. 
    {
        return triangleToCell_; 
    }

    const auto& vertexToCell() const 
    {
        return vertexToCell_; 
    }

    //- Interface cell -> contained triangles
    const auto& interfaceCellToTriangles() const
    {
        return interfaceCellToTriangles_;
    }

    //- Interface cell -> contained vertices
    const auto& interfaceCellToVertices() const
    {
        return interfaceCellToVertices_;
    }

    template<typename LabelList> 
    void setTriangleToCell(LabelList const& l)
    {
        triangleToCell_ = l; 
    }

    //- Cell -> Nearest Triangle.  
    const DynamicList<pointIndexHit>& cellsTriangleNearest() const
    {
        return cellsTriangleNearest_; 
    };

    DynamicList<pointIndexHit>& cellsTriangleNearest()
    {
        return cellsTriangleNearest_; 
    };

    //- Point -> Nearest Triangle.  
    DynamicList<pointIndexHit>& pointsTriangleNearest()
    {
        return pointsTriangleNearest_; 
    };

    const DynamicList<pointIndexHit>& pointsTriangleNearest() const
    {
        return pointsTriangleNearest_; 
    };


    // Update
    
    void update();  

    void updateVertexToCell();  

    void updateInterfaceCellToTriangles();
    void updateInterfaceCellToVertices();

    // Empty writer: regIOobject is used only for registry queries.
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
