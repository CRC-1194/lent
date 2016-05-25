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


#include "lentCommunication.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(lentCommunication, 0);
    defineRunTimeSelectionTable(lentCommunication, FrontMesh);
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
    const word name = configDict.lookup("type");

    FrontMeshConstructorTable::iterator cstrIter =
        FrontMeshConstructorTablePtr_->find(name);

    if (cstrIter == FrontMeshConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "lentCommunication::New(const word& name)"
        )   << "Unknown lentCommunication type "
            << name << nl << nl
            << "Valid lentCommunications are : " << endl
            << FrontMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<lentCommunication> (cstrIter()(front, mesh));
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
                if (foundCell > 0)
                {
                    // Set the triangle cell to the found cell.
                    triangleToCell_[triangleI] = foundCell;
                    // Set the vertex cell to the found cell.
                    vertexToCell_[triangle[vertexI]] = foundCell;
                }
            }
        }
    }
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
}

bool lentCommunication::writeData(Ostream& os) const
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
