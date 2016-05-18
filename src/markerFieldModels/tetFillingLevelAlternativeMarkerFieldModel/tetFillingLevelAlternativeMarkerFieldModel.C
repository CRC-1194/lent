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
    Foam::tetFillingLevelAlternativeMarkerFieldModel

SourceFiles
    tetFillingLevelAlternativeMarkerFieldModel.C

Author
    Tobias Tolle    tolle@csi.tu-darmstadt.de

Description
    Computes the markerfield from the vertex based and cell centre based 
    signed distances by decomposition of a cell into tetrahedrons and calculating
    the filling level for each tetrahedron as described in

        "From level set to volume of fluid and back again at second-order 
         accuracy" Detrixhe and Aslam

    This class uses an alternative approach for the triangulation of the faces
    which does not introduce additional points and thus avoids interpolation
    errors regardign the signed distance.

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

#include "addToRunTimeSelectionTable.H"

#include "tetFillingLevelAlternativeMarkerFieldModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(tetFillingLevelAlternativeMarkerFieldModel, 0);
    addToRunTimeSelectionTable(markerFieldModel, tetFillingLevelAlternativeMarkerFieldModel, Dictionary);

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
List<tetrahedron> tetFillingLevelAlternativeMarkerFieldModel::tetDecomposition(const point& cellCentre,
      const face& cellFace,
      const pointField& points,
      const pointScalarField& pointDistance,
      scalar centreDistance) const
{
    // Perform decomposition by using the first point of the face and
    // combining it with with each edge which is no direct neighbour of
    // this point. The result are nEdges-2 triangles.
    List<tetrahedron> tets(cellFace.nEdges()-2);

    forAll(tets, I)
    {
        tets[I].distance.resize(4);
    }

    edgeList edges = cellFace.edges();

    for(label edgeI=1; edgeI < (edges.size()-1); edgeI++)
    {
        edge e = edges[edgeI];

        tets[edgeI-1].vertex[0] = cellCentre;
        tets[edgeI-1].vertex[1] = points[cellFace[0]];
        tets[edgeI-1].vertex[2] = points[e[0]];
        tets[edgeI-1].vertex[3] = points[e[1]];

        tets[edgeI-1].distance[0] = centreDistance;
        tets[edgeI-1].distance[1] = pointDistance[cellFace[0]];
        tets[edgeI-1].distance[2] = pointDistance[e[0]];
        tets[edgeI-1].distance[3] = pointDistance[e[1]];
    }

    return tets;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
