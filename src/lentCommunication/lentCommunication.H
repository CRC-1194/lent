/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
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

Class
    Foam::lentCommunication

SourceFiles
    lentCommunication.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
        Front / Mesh communication maps.  

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef lentCommunication_H
#define lentCommunication_H

#include "dictionary.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "DynamicList.H"
#include "frontMeshSearch.H"
#include "triSurfaceFront.H"
#include "fvMesh.H"

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
    );

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
    const auto& triangleToCell() const 
    {
        return triangleToCell_; 
    }

    auto& triangleToCell() // TODO: Remove after removing mesh search from the velocity interpolator. Follow through in lentMethod class. TM. 
    {
        return triangleToCell_; 
    }

    const auto& vertexToCell() const 
    {
        return vertexToCell_; 
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
