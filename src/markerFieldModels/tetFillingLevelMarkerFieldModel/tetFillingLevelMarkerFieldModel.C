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
    Foam::tetFillingLevelMarkerFieldModel

SourceFiles
    tetFillingLevelMarkerFieldModel.C

Author
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Computes the markerfield from the vertex based and cell centre based 
    signed distances by decomposition of a cell into tetrahedrons and calculating
    the filling level for each tetrahedron as described in

        "From level set to volume of fluid and back again at second-order 
         accuracy" Detrixhe and Aslam

    This class uses a barycentric decomposition of the cell faces into triangles.

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

#include <utility>
#include <cassert>

#include "addToRunTimeSelectionTable.H"
#include "fvcAverage.H"
#include "surfaceInterpolate.H"
#include "volFields.H"

#include "tetFillingLevelMarkerFieldModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(tetFillingLevelMarkerFieldModel, 0);
    addToRunTimeSelectionTable(markerFieldModel, tetFillingLevelMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //
tetFillingLevelMarkerFieldModel::tetFillingLevelMarkerFieldModel(
                                    const dictionary& configDict)
:
    markerFieldModel{configDict},
    pointDistFieldName_{configDict.get<word>("pointDistance")},
    nSmoothingSteps_{configDict.get<label>("nSmoothingSteps")}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void tetFillingLevelMarkerFieldModel::calcMarkerField(volScalarField& markerField) const
{
    const fvMesh& mesh = markerField.mesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbours = mesh.faceNeighbour();
    const faceList& faces = mesh.faces();

    tagInterfaceCells(markerField);

    label cellI = -1;

    forAll(faces, faceI)
    {
        // Check if owner or neighbour cell or both are an interface cell.
        // If true, perform the barycentric decomposition and computation
        // of filling level contribution according to the paper mentioned in
        // the class description
        if (markerField[owners[faceI]] >= 0.0)
        {
            cellI = owners[faceI];
            markerField[cellI] += fillingLevel(cellI, faces[faceI], mesh);
        }

        if (faceI < neighbours.size())
        {
            if (markerField[neighbours[faceI]] >= 0.0)
            {
               cellI = neighbours[faceI];
               markerField[cellI] += fillingLevel(cellI, faces[faceI], mesh);
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
    smoothMarkerField(markerField);
    cutOverUnderShoots(markerField);
}

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
        // Check for LARGE to ensure only actual interface cells are tagged
        // --> makes this approach more resilient against errors in the 
        // signed distance propagation hopefully
        if ((firstPointDist * secondPointDist ) < 0
            && mag(firstPointDist) < GREAT 
            && mag(secondPointDist) < GREAT
           )
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
    }
}

List<tetrahedron> tetFillingLevelMarkerFieldModel::tetDecomposition
(
    const point& cellCentre,
    const face& cellFace,
    const pointField& points,
    const pointScalarField& pointDistance,
    scalar centreDistance
) const
{
    // Perform decomposition using the barycentre of the face
    List<tetrahedron> tets(cellFace.nEdges());

    forAll(tets, I)
    {
        tets[I].distance.resize(4);
    }

    point faceCentre = cellFace.centre(points);
    scalar faceDistance = cellFace.average(points, pointDistance);
    edgeList edges = cellFace.edges();

    forAll(edges, edgeI)
    {
        edge e = edges[edgeI];

        tets[edgeI].vertex[0] = cellCentre;
        tets[edgeI].vertex[1] = faceCentre;
        tets[edgeI].vertex[2] = points[e[0]];
        tets[edgeI].vertex[3] = points[e[1]];

        tets[edgeI].distance[0] = centreDistance;
        tets[edgeI].distance[1] = faceDistance;
        tets[edgeI].distance[2] = pointDistance[e[0]];
        tets[edgeI].distance[3] = pointDistance[e[1]];
    }

    return tets;
}

scalar tetFillingLevelMarkerFieldModel::fillingLevel
(
    const label& cellID, const face& cellFace, const fvMesh& mesh
) const
{
    scalar absoluteFilling = 0.0;

    const volScalarField& signedDistance =
        mesh.lookupObject<volScalarField>(markerFieldModel::cellDistFieldName());
    const pointScalarField& pointDistance =
        mesh.lookupObject<pointScalarField>(pointDistFieldName()); 
    const pointField& points = mesh.points();
    const edgeList edges = cellFace.edges();

    point cellCentre = mesh.C()[cellID];

    List<tetrahedron> tets = tetDecomposition(cellCentre, cellFace, points, pointDistance, signedDistance[cellID]);

    forAll(tets, I)
    {
        absoluteFilling +=
            tetrahedralVolume(tets[I]) * volumeFraction(tets[I].distance);
    }

    assert(absoluteFilling < mesh.V()[cellID]);
    return absoluteFilling;
}

scalar tetFillingLevelMarkerFieldModel::tetrahedralVolume(const tetrahedron& tet) const
{
    return fabs((tet.vertex[1] - tet.vertex[0]) & 
           ((tet.vertex[2] - tet.vertex[0]) ^ (tet.vertex[3] - tet.vertex[0])))
           /6.0;
}

scalar tetFillingLevelMarkerFieldModel::volumeFraction(SortableList<scalar>& d) const
{
    // This function implements the computation of the volume fraction from
    // the signed distance for a tetrahedron as described by the paper
    // mentioned in the class description
    assert(d.size() == 4);
    
    scalar volFraction = 0.0;
    d.sort();
    label negativeEntries = countNegativeEntries(d);

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
        volFraction = ( d[0]*d[1] * (d[2]*d[2] + d[2]*d[3] + d[3]*d[3])
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

    assert (volFraction >= 0.0 && volFraction <= 1.0);
    return volFraction;
}

label tetFillingLevelMarkerFieldModel::countNegativeEntries(List<scalar>& distance) const
{
    label negativeEntries = 0;

    forAll(distance, I)
    {
        if (distance[I] <= 0) negativeEntries++;
    }

    return negativeEntries;
}

void tetFillingLevelMarkerFieldModel::smoothMarkerField(volScalarField& markerField) const
{
    for (label I = 0; I < nSmoothingSteps_; I++)
    {
        markerField = fvc::average(fvc::interpolate(markerField));
    }
}

void tetFillingLevelMarkerFieldModel::cutOverUnderShoots(volScalarField& markerField) const
{
    forAll(markerField, cellI)
    {
        if (markerField[cellI] < 0.0)
        {
            Info << "Alpha undershoot in cell " << cellI
                 << ", value = " << markerField[cellI] << endl;
            markerField[cellI] = 0.0;
        }
        else if (markerField[cellI] > 1.0)
        {
            Info << "Alpha overshoot in cell " << cellI
                 << ", value = " << markerField[cellI] << endl;
            markerField[cellI] = 1.0;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
