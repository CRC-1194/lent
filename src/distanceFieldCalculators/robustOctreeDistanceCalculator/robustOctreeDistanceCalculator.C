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
    Foam::robustOctreeDistanceCalculator

SourceFiles
    robustOctreeDistanceCalculator.H

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Computes signed distance fields between the volume mesh and the immersed
    surface mesh using octree searches implemented in the optimizedOctree class.

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


#include "robustOctreeDistanceCalculator.H"
#include "addToRunTimeSelectionTable.H"
#include "volumeType.H"
#include "lentCommunication.H"

// From SMCIA
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "triSurfaceTools.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(robustOctreeDistanceCalculator, 0);
    addToRunTimeSelectionTable(
       distanceFieldCalculator,
       robustOctreeDistanceCalculator,
       Dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
vectorField robustOctreeDistanceCalculator::computeVertexNormals(const triSurface& front) const
{
    vectorField vertexNormals{front.nPoints(), vector{0,0,0}};

    /* The vertex normals are computed as a weighted sum of normals
     * of the adjacent triangles. The weights are the triangle angles at the
     * vertex at hand.
     * See
     *      "Computing Vertex Normals from Polygonal Facets",
     *      G. Thürmer & C. Wüthrich (2012)
     *      https://doi.org/10.1080/10867651.1998.10487487
     * for details. There is a proof that this gives correct inside/outside
     * information by J. Baerentzen and H. Aanaes, Technical University of
     * Denmark
     */
    const auto& tri_normals = front.faceNormals();
    const auto& vertex_to_faces = front.pointFaces();
    const auto& vertices = front.localPoints();

    forAll(vertices, v_id)
    {
        label vid_a{v_id};
        label vid_b{0};
        label vid_c{0};

        for (const auto fid : vertex_to_faces[v_id])
        {
            triSurfaceTools::otherVertices(front, fid, vid_a, vid_b, vid_c);
            vector v1{vertices[vid_b] - vertices[vid_a]};
            vector v2{vertices[vid_c] - vertices[vid_a]};
            scalar alpha{Foam::acos(
                std::clamp((v1 & v2) / (mag(v1) * mag(v2)), -1.0, 1.0))};

            vertexNormals[v_id] += alpha * tri_normals[fid];
        }

        vertexNormals[v_id] /= mag(vertexNormals[v_id]) + SMALL;
    }

    return vertexNormals;
}
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

robustOctreeDistanceCalculator::robustOctreeDistanceCalculator(
    const dictionary& config
)
:
    distanceFieldCalculator(config),
    narrowBandTmp_(
       narrowBandPropagation::New(config.subDict("narrowBandPropagation"))
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void robustOctreeDistanceCalculator::calcCellsToFrontDistance(
    volScalarField& signedDistance,
    const volScalarField& searchDistanceSqr,
    const triSurfaceFront& front
)
{
    signedDistance = dimensionedScalar("GREAT", dimLength, GREAT);

    const auto& vertexNormals = computeVertexNormals(front);

    // Get the lent communication structure from the registry. 
    // Get the non-const reference because the nearest maps are modified by 
    // the distance calculation. TM. 
    const fvMesh& mesh = signedDistance.mesh();
    lentCommunication& comm = const_cast<lentCommunication&>(
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
        )
    );

    const volVectorField& C = mesh.C();
    auto& cellsTriangleNearest = comm.cellsTriangleNearest(); 

    triSurfaceSearch surfaceSearch{front};
    surfaceSearch.findNearest(C, searchDistanceSqr, cellsTriangleNearest);

    forAll(cellsTriangleNearest, I)
    {
        const pointIndexHit& hitInfo = cellsTriangleNearest[I];

        if (hitInfo.hit())
        {
            vector delta_v{C[I] - hitInfo.hitPoint()};
            auto snormal = normalAtSurface(hitInfo, front, vertexNormals);
            signedDistance[I] = mag(delta_v) * sign(snormal & delta_v);
        }
    }

    narrowBandTmp_->ensureNarrowBand(signedDistance, GREAT);
}

void robustOctreeDistanceCalculator::calcPointsToFrontDistance(
    pointScalarField& pointSignedDistance,
    const pointScalarField& pointSearchDistanceSqr,
    const triSurfaceFront& front
)
{
    pointSignedDistance = dimensionedScalar("GREAT", dimLength, GREAT);

    const auto& vertexNormals = computeVertexNormals(front);

    // Get the lent communication structure from the registry. 
    // Get the non-const reference because the nearest maps are modified by 
    // the distance calculation. TM. 
    const polyMesh& mesh = pointSignedDistance.mesh()();
    lentCommunication& comm = const_cast<lentCommunication&>(
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
        )
    );

    // Get the cell centres.
    const pointMesh& pMesh = pointSignedDistance.mesh();
    const pointField& points = pMesh().points();
    auto& pointsTriangleNearest = comm.pointsTriangleNearest(); 

    triSurfaceSearch surfaceSearch{front};
    surfaceSearch.findNearest(points, pointSearchDistanceSqr, pointsTriangleNearest);

    forAll(pointsTriangleNearest, I)
    {
        const pointIndexHit& hitInfo = pointsTriangleNearest[I];

        if (hitInfo.hit())
        {
            vector delta_v{points[I] - hitInfo.hitPoint()};
            auto snormal = normalAtSurface(hitInfo, front, vertexNormals);
            pointSignedDistance[I] = mag(delta_v) * sign(snormal & delta_v);
        }
    }

    // TODO: GREAT --> make it a controllable configuraton value
    narrowBandTmp_->ensureNarrowBand(pointSignedDistance, GREAT);
}

vector robustOctreeDistanceCalculator::normalAtSurface(
    const pointIndexHit& hitInfo,
    const triSurface front,
    const vectorField& vertexNormals
) const
{
    vector normal{0, 0, 0};

    const auto& fnormals = front.faceNormals();
    const auto& v = front.localPoints();

    // Transformation to local coordinate system:
    // - origin: point a (first point of triangle)
    // - x-axis: point a to point b (second point of triangle)
    // - y-axis: cross product of triangle normal and ex
    const auto tri_hit = front.localFaces()[hitInfo.index()];
    const auto a_to_b = v[tri_hit[1]] - v[tri_hit[0]];
    const auto a_to_c = v[tri_hit[2]] - v[tri_hit[0]];
    const auto a_to_hit = hitInfo.hitPoint() - v[tri_hit[0]];
    const auto ex = a_to_b / mag(a_to_b);
    const auto ey = fnormals[hitInfo.index()] ^ ex;
    const vector b{a_to_b & ex, 0, 0};
    const vector c{a_to_c & ex, a_to_c & ey, 0};
    const vector h_l{a_to_hit & ex, a_to_hit & ey, 0};

    // This is a linear shape function for a triangle where the vertices are:
    // - 1) the origin (0,0)
    // - 2) a point on the x-axis (t,0)
    // - 3) an arbitrary point (u,v)
    normal = vertexNormals[tri_hit[0]] *
            (1.0 - h_l.x() / b.x() + h_l.y() / c.y() * (c.x() / b.x() - 1.0)) +
        vertexNormals[tri_hit[1]] *
            (h_l.x() / b.x() - h_l.y() * c.x() / (b.x() * c.y())) +
        vertexNormals[tri_hit[2]] * (h_l.y() / c.y());

    return normal;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
