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


#include "lentCommunication.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(lentCommunication, 0);
    defineRunTimeSelectionTable(lentCommunication, FrontMesh)
    addToRunTimeSelectionTable(lentCommunication, lentCommunication, FrontMesh);

    word lentCommunication::registeredName(
            const triSurfaceFront& front, 
            const polyMesh& mesh
    )
    {
        return mesh.name() + "-" + front.name() + "-communication"; 
    } 


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentCommunication::lentCommunication(
    const triSurfaceFront& front, 
    const fvMesh& mesh
)
    :
        regIOobject(
            IOobject(
               lentCommunication::registeredName(front,mesh), 
               mesh.thisDb().instance(),  
               mesh,
               IOobject::NO_READ, 
               IOobject::NO_WRITE
            )
        ),
        front_(front), 
        mesh_(mesh), 
        searchAlg_(),
        triangleToCell_(front_.nFaces()),
        vertexToCell_(front_.nPoints()),
        interfaceCellToTriangles_{},
        interfaceCellToVertices_{},
        cellsTriangleNearest_(), 
        pointsTriangleNearest_()
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<lentCommunication>
lentCommunication::New(
        const dictionary& configDict, 
        const triSurfaceFront& front, 
        const fvMesh& mesh
)
{
    const word name = configDict.get<word>("type");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = FrontMeshConstructorTable(name);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "lentCommunication",
            name,
            *FrontMeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object.
    return autoPtr<lentCommunication>(ctorPtr(front, mesh));
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Updates both triangle->cell and vertex->cell maps using the KVS search algorithm. 
// Used after front evolution, does not account for topological changes of the front.
void lentCommunication::update()
{
    const auto& triangles = front_.localFaces();
    const auto& vertices = front_.points();

    forAll (triangleToCell_, triangleI) 
    {
        const auto& triangle = triangles[triangleI];

        forAll (triangle, vertexI)
        {
            label foundCell = -1;

            const point& vertex = vertices[triangle[vertexI]]; 

            // If the vertex is within a the triangleToCell cell. 
            if (searchAlg_.pointIsInCell(vertex, triangleToCell_[triangleI], mesh_)) 
            {
                // Existing cell found.
                foundCell = triangleToCell_[triangleI];

                // DEBUGGING
                //WarningInFunction 
                    //<< "Existing cell " 
                    //<< triangleToCell_[triangleI] 
                    //<< " and front vertex " 
                    //<< triangle[vertexI] << endl; 
                    
                // Set the vertex cell to the same cell.
                vertexToCell_[triangle[vertexI]] = triangleToCell_[triangleI];
            } else
            {
                // Find the cell that contains the vertex.
                foundCell = searchAlg_.cellContainingPoint(
                    vertex,
                    mesh_,
                    triangleToCell_[triangleI]
                );

                // If the cell is found. 
                if (foundCell != -1)
                {
                    // Set the triangle cell to the found cell.
                    triangleToCell_[triangleI] = foundCell;
                    // Set the vertex cell to the found cell.
                    vertexToCell_[triangle[vertexI]] = foundCell;

                    // DEBUGGING
                    //WarningInFunction 
                        //<< "KVS found cell " 
                        //<< foundCell  
                        //<< " and front vertex " 
                        //<< triangle[vertexI] << endl; 
                }
            }

            // DEBUGGING
            //if (foundCell == -1)
                //WarningInFunction 
                    //<< "Front vertex not found for cell " 
                    //<< triangleToCell_[triangleI] 
                    //<< " and front vertex " 
                    //<< triangle[vertexI] << endl; 
        }
    }

    updateInterfaceCellToTriangles();
    updateInterfaceCellToVertices();
}

// Reconstruction results in a triangle->cell relationship, regardless which  
// reconstruction algorithm is used.  
// Use triangle->cell to update vertex->cell map.
void lentCommunication::updateVertexToCell()
{
    const auto& triangles = front_.localFaces();
    const auto& vertices = front_.points();

    vertexToCell_.resize(vertices.size(), -1); 

    forAll (triangleToCell_, triangleI) 
    {
        const auto& triangle = triangles[triangleI];

        forAll (triangle, vertexI)
        {
            const point& vertex = vertices[triangle[vertexI]]; 
            // If the vertex is within a the triangleToCell cell. 
            if (searchAlg_.pointIsInCell(vertex, triangleToCell_[triangleI], mesh_)) 
            {
                // Set the vertex cell to the found cell.
                vertexToCell_[triangle[vertexI]] = triangleToCell_[triangleI];
            }
            else
            {
                // Find the cell that contains the vertex.
                const label foundCell = searchAlg_.cellContainingPoint(
                    vertex,
                    mesh_,
                    triangleToCell_[triangleI]
                );

                // If the cell is found. 
                if (foundCell > 0)
                {
                    // Set the vertex cell to the found cell.
                    vertexToCell_[triangle[vertexI]] = foundCell;
                }
            }
        }
    }

    updateInterfaceCellToVertices();
}

void lentCommunication::updateInterfaceCellToTriangles()
{
    interfaceCellToTriangles_.clear();

    forAll(triangleToCell_, I)
    {
        interfaceCellToTriangles_[triangleToCell_[I]].push_back(I);
    }
}

void lentCommunication::updateInterfaceCellToVertices()
{
    interfaceCellToVertices_.clear();

    forAll(vertexToCell_, I)
    {
        interfaceCellToVertices_[vertexToCell_[I]].push_back(I);
    }
}

bool lentCommunication::writeData(Ostream&) const
{
    FatalErrorIn("lentMethod::writeData(Ostream& os)")
    << "lentMethod is not supposed to be written "
    << "regIOobject inherited to allow registry queries." << endl;

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
