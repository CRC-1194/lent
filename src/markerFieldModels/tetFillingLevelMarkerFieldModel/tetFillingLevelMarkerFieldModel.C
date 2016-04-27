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
    Foam::tetFillingLevelMarkerFieldModel

SourceFiles
    tetFillingLevelMarkerFieldModel.C

Author
    Tobias Tolle    tolle@csi.tu-darmstadt.de

Description
    Computes the markerfield from the vertex based and cell centre based 
    signed distances by decomposition of a cell into tetrahedrons and calculating
    the filling level for each tetrahedron as described in

        "From level set to volume of fluid and back again at second-order 
         accuracy" Detrixhe and Aslam

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

#include <utility>
#include <cassert>

#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

#include "tetFillingLevelMarkerFieldModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(tetFillingLevelMarkerFieldModel, 0);
    addToRunTimeSelectionTable(markerFieldModel, tetFillingLevelMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void tetFillingLevelMarkerFieldModel::tagInterfaceCells(volScalarField& markerField) const
{
    // For clear distinction, set marker field value of bulk cells
    // temporarily to a negative value
    forAll(markerField, cellI)
    {
        markerField[cellI] = -GREAT;
    }

    const fvMesh& mesh = markerField.mesh(); 
    const edgeList& meshEdges = mesh.edges(); 
    const labelListList& meshEdgeCells = mesh.edgeCells();  

    const pointScalarField& pointDistance = 
        mesh.lookupObject<pointScalarField>(pointDistFieldName()); 

    // Find an edge of a cell where the point signed distance switches in sign.
    // Find the cells that share that edge and set the marker field to zero
    // --> tagged as interface cell and ready for further calculations
    forAll(meshEdges, edgeI)
    {
        const edge& meshEdge = meshEdges[edgeI]; 
        auto firstPointDist = pointDistance[meshEdge[0]];
        auto secondPointDist = pointDistance[meshEdge[1]];

        // If the distance switches sign for an edge of the cell cellI. 
        if ((firstPointDist * secondPointDist ) < 0)
        {
            const labelList& edgeCells = meshEdgeCells[edgeI]; 
            forAll(edgeCells, edgeJ)
            {
                const label cellLabel = edgeCells[edgeJ]; 

                markerField[cellLabel] = 0.0;
            }
        }
    }
}

void tetFillingLevelMarkerFieldModel::setBulkMarkerField(volScalarField& markerField) const
{
    const fvMesh& mesh = markerField.mesh();

    const volScalarField& signedDistance =
        mesh.lookupObject<volScalarField>(markerFieldModel::cellDistFieldName());

    forAll(markerField, cellI)
    {
        if (markerField[cellI] == -GREAT)
        {
            if (signedDistance[cellI] < 0.0)
            {
                markerField[cellI] = 0.0;
            }
            else 
            {
                markerField[cellI] = 1.0;
            }
        }

        assert (markerField[cellI] >= 0.0);
    }
}

scalar tetFillingLevelMarkerFieldModel::fillingLevel
(
    const label& cellID, const face& cellFace, const fvMesh& mesh
) const
{
    const volScalarField& signedDistance =
        mesh.lookupObject<volScalarField>(markerFieldModel::cellDistFieldName());
    const pointScalarField& pointDistance =
        mesh.lookupObject<pointScalarField>(pointDistFieldName()); 
    const pointField& points = mesh.points();
    const edgeList edges = cellFace.edges();

    point cellCentre = mesh.C()[cellID];
    point faceCentre = cellFace.centre(points);

    // Stores the signed distances for a tetrahedron
    List<scalar> distances(4);

    scalar absoluteFilling = 0.0;
    scalar distanceFaceCentre = 0.0;

    // For now, compute the signed distance at the face centre by the
    // arithmetic mean of the vertex values
    // TODO: check quality and replace with more accurate scheme if required (TT)
    forAll(cellFace, vertexI)
    {
        distanceFaceCentre += pointDistance[cellFace[vertexI]];
    }
    distanceFaceCentre /= cellFace.size();

    forAll(edges, edgeI)
    {
        distances[0] = signedDistance[cellID];
        distances[1] = distanceFaceCentre;
        distances[2] = pointDistance[edges[edgeI][0]];
        distances[3] = pointDistance[edges[edgeI][1]];

        absoluteFilling +=
            tetrahedralVolume
            (
                cellCentre, faceCentre,
                points[edges[edgeI][0]], points[edges[edgeI][1]]
            ) 
            * volumeFraction(distances);
    }

    return absoluteFilling;
}

scalar tetFillingLevelMarkerFieldModel::tetrahedralVolume
(
    const point& pointA, const point& pointB, const point& pointC,
    const point& pointD
) const
{
    return fabs((pointA - pointB) & ((pointC - pointB) ^ (pointD - pointB)))/6.0;
}


scalar tetFillingLevelMarkerFieldModel::volumeFraction(List<scalar>& d) const
{
    // This function implements the computation of the volume fraction from
    // the signed distance for a tetrahedron as described by the paper
    // mentioned in the class description
    scalar volFraction = 0.0;
    label negativeEntries = sortDistances(d);

    if (negativeEntries == 4)
    {
        volFraction = 0.0;
    }
    else if (negativeEntries == 3)
    {
        volFraction = std::pow(d[3],3) /
               ((d[3] - d[0]) * (d[3] - d[1]) * (d[3] - d[2]));
    }
    else if (negativeEntries == 2)
    {
        volFraction = 1.0 - ( d[0]*d[1] * (d[2]*d[2] + d[2]*d[3] + d[3]*d[3])
                + d[2]*d[3] * (d[2]*d[3] - (d[0]+d[1])*(d[2]+d[3])) )
                / ((d[0]-d[2]) * (d[1]-d[2]) * (d[0]-d[3]) * (d[1]-d[3]));
    }
    else if (negativeEntries == 1)
    {
        volFraction = 1.0 + std::pow(d[0],3) / 
                ( (d[1]-d[0]) * (d[2]-d[0]) * (d[3]-d[0]) );
    }
    else
    {
        volFraction = 1.0;
    }

    assert (volFraction >= 0.0);
    return volFraction;
}

label tetFillingLevelMarkerFieldModel::sortDistances(List<scalar>& distances) const
{
    label negativeEntries = 0;

    forAll(distances, I)
    {
        for (label K = I+1; K < 4; K++)
        {
            if (distances[K] < distances[I])
            {
                std::swap(distances[I], distances[K]);
            }
        }

        if (distances[I] <= 0.0)
        {
            negativeEntries++;
        }
    }

    return negativeEntries;
}

// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
tetFillingLevelMarkerFieldModel::tetFillingLevelMarkerFieldModel(
                                    const dictionary& configDict)
:
    markerFieldModel(configDict),
    pointDistFieldName_(configDict.lookup("pointDistance"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void tetFillingLevelMarkerFieldModel::calcMarkerField(volScalarField& markerField) const
{
    tagInterfaceCells(markerField);

    const fvMesh& mesh = markerField.mesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbours = mesh.faceNeighbour();
    const faceList& faces = mesh.faces();

    label cellI = -1;

    forAll(faces, faceI)
    {
        // Check if either owner or neighbour cell is an interface cell
        // If true, perform the barycentric decomposition and computation
        // of filling level contribution according to the paper mentioned in
        // the class description
        if (markerField[owners[faceI]] >= 0.0)
        {
            cellI = owners[faceI];
            markerField[cellI] += fillingLevel(cellI, faces[faceI], mesh);
            cellI = -1;
        }

        if (faceI < neighbours.size())
        {
            if (markerField[neighbours[faceI]] >= 0.0)
            {
               cellI = neighbours[faceI];
               markerField[cellI] += fillingLevel(cellI, faces[faceI], mesh);
               cellI = -1;
            }
        }
    }

    // Normalize value of interface cells to acquire actual volume fractions
    forAll(markerField, cellI)
    {
        if (markerField[cellI] > 0.0)
        {
            markerField[cellI] /= mesh.V()[cellI];
        }
    }

    setBulkMarkerField(markerField);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
