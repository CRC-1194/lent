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
    Foam::FrontTracking::frontSmoother

SourceFiles
    frontSmoother.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description
    Implementation of two volume-conservative smoothing algorithms for
    triangulated surfaces as proposed by Kuprat et al 2001
    Both algorithms have been extended so they are applicable to fronts
    with boundaries, e.g. in 2D cases or contact line problems

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

#include <algorithm>

#include "lentCommunication.H"

#include "frontSmoother.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void frontSmoother::updateGlobalPoints(triSurfaceFront& front) const
{
    pointField& frontPoints = const_cast<pointField&>(front.points());
    const auto& localPoints = front.localPoints();
    const auto& globalToLocal = front.meshPoints();

    forAll(globalToLocal, pointI)
    {
        frontPoints[globalToLocal[pointI]] = localPoints[pointI];
    }
}

vector frontSmoother::computeA(const label& pointLabel, const triSurfaceFront& front) const
{
    vector A{0,0,0};

    const auto& connectedFaces = front.pointFaces()[pointLabel];
    const auto& faces = front.localFaces();
    const auto& points = front.localPoints();

    forAll(connectedFaces, index)
    {
        const auto& aFace = faces[connectedFaces[index]];
        A += (points[aFace[1]] - points[aFace[0]])
                    ^ (points[aFace[2]] - points[aFace[0]]);
    }

    return A;
}

vector frontSmoother::computeV(const label& edgeLabel, const triSurfaceFront& front) const
{
    vector V{0,0,0};

    const auto& points = front.localPoints();
    const auto& relaxEdge = front.edges()[edgeLabel];
    const auto& connectedFaces = front.edgeFaces()[edgeLabel];
    const auto& faces = front.localFaces();

    forAll(connectedFaces, index)
    {
        V *= -1.0;

        const auto& aFace = faces[connectedFaces[index]];

        forAll(aFace, I)
        {
            if (aFace[I] != relaxEdge[0] && aFace[I] != relaxEdge[1])
            {
                V += points[aFace[I]];
            }
        }
    }

    // Check orientation
    const auto& x0 = points[relaxEdge[0]];
    const auto& x1 = points[relaxEdge[1]];
    const auto& aFace = faces[connectedFaces[0]];
    if ((aFace.normal(points) & ((x0 - x1) ^ V)) < 0.0)
    {
        V *= -1.0;
    }

    return V;
}

vectorTuple frontSmoother::smoothEdge(const edge& relaxEdge, const triSurfaceFront& front) const
{
    vectorTuple smoothedPoints{vector{0,0,0}, vector{0,0,0}};

    smoothedPoints[0] = starBarycentre(relaxEdge[0], front);
    smoothedPoints[1] = starBarycentre(relaxEdge[1], front);

    return smoothedPoints;
}

vector frontSmoother::starBarycentre(const label& pointLabel, const triSurfaceFront& front) const
{
    vector starBC{0,0,0};
    scalar starArea{0};

    const auto& connectedFaces = front.pointFaces()[pointLabel];
    const auto& faces = front.localFaces();
    const auto& points = front.localPoints();

    forAll(connectedFaces, index)
    {
        const auto& aFace = faces[connectedFaces[index]];

        starBC += aFace.centre(points) * aFace.mag(points);
        starArea += aFace.mag(points);
    }

    starBC /= starArea;

    return starBC;
}

std::vector<label> frontSmoother::singlePointBoundaryEdges(const triSurfaceFront& front) const
{
    // Remember: OpenFOAM only counts edges with both points located
    // on the boundary as boundary edges
    std::vector<label> singlePointBoundaryEdges{};

    const auto& boundaryPoints = front.boundaryPoints();
    const auto& pointToEdge = front.pointEdges();
    auto nInternalEdges = front.nInternalEdges();

    forAll(boundaryPoints, index)
    {
        auto connectedEdges = pointToEdge[boundaryPoints[index]];

        forAll(connectedEdges, K)
        {
            if (connectedEdges[K] < nInternalEdges)
            {
                singlePointBoundaryEdges.push_back(connectedEdges[K]);
            }
        }
    }

    std::sort(singlePointBoundaryEdges.begin(), singlePointBoundaryEdges.end());

    return singlePointBoundaryEdges;
}

std::vector<label> frontSmoother::internalEdges(const label nInternalEdges, const std::vector<label>& boundaryEdges) const
{
    std::vector<label> internalEdges{};

    label boundaryEdgeIndex = 0;

    for (label index = 0; index < nInternalEdges; ++index)
    {
        if (!boundaryEdges.empty() && index == boundaryEdges[boundaryEdgeIndex])
        {
            boundaryEdgeIndex++;
            continue;
        }
        else
        {
            internalEdges.push_back(index);
        }
    }

    return internalEdges;
}

label frontSmoother::determineBoundaryPoint(const edge& relaxEdge, const triSurfaceFront& front) const
{
    const auto& boundaryPoints = front.boundaryPoints();

    if (std::binary_search(boundaryPoints.begin(), boundaryPoints.end(), relaxEdge[0]))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

Tensor<scalar> frontSmoother::movementRestriction(const label& pointLabel, const triSurfaceFront& front, const fvMesh& mesh) const
{
    Tensor<scalar> restriction{};

    // Find face on which the point is located
    auto minDistFaceLabel = containingFace(pointLabel, front, mesh);

    // Check if the point coincides with an edge
    auto edgeWithFrontVertex = containingEdge(pointLabel, minDistFaceLabel, front, mesh);

    // If front vertex coincides with an edge: are both connected faces coplanar?
    bool facesAreCoplanar = false;

    if (edgeWithFrontVertex > -1)
    {
        facesAreCoplanar = boundaryFacesAreCoplanar(edgeWithFrontVertex, front, mesh);
    }

    if (edgeWithFrontVertex < 0 || facesAreCoplanar)
    {
        auto faceNormal = mesh.faceAreas()[minDistFaceLabel];
        faceNormal /= mag(faceNormal);
        restriction = Identity<scalar>{} - faceNormal*faceNormal;
    }
    else
    {
        const auto& frontVertexEdge = front.edges()[edgeWithFrontVertex];
        auto edgeDirection = frontVertexEdge.vec(front.localPoints());
        edgeDirection /= mag(edgeDirection);
        restriction = edgeDirection*edgeDirection;
    }

    return restriction;
}

label frontSmoother::containingFace(const label& pointLabel, const triSurfaceFront& front, const fvMesh& mesh) const
{
    label minDistFaceLabel = -1;

    const auto& frontVertex = front.localPoints()[pointLabel];
    const auto& globalPointLabel = front.meshPoints()[pointLabel];
    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(front,mesh)
                                ); 
    const auto& containingCellLabel = communication.vertexToCell()[globalPointLabel];
    const auto& containingCell = mesh.cells()[containingCellLabel];
    const auto& faceCentre = mesh.faceCentres();
    auto nInternalFaces = mesh.nInternalFaces();

    // Find closest face centre to point --> thats the boundary face we are
    // looking for
    scalar minDist{GREAT};

    forAll(containingCell, index)
    {
        // Only boundary faces are viable for computing a boundary normal
        if (containingCell[index] < nInternalFaces)
        {
            continue;
        }
        // Assumption: a boundary point of the front is always located at
        // the boundary of the fvMesh, either because of a 2D simulation
        // or there is a contact line. Therefore, this boundary point must
        // be located in a plane spanned by a boundary face. Thus, the distance
        // between the boundary point and the corresponding face centre should
        // be zero if projected to the face normal
        auto vToFaceCentre = frontVertex - faceCentre[containingCell[index]];
        auto projectedDistance = mag(vToFaceCentre & mesh.faceAreas()[containingCell[index]] / mag(mesh.faceAreas()[containingCell[index]]));

        if (minDist > projectedDistance)
        {
            minDist = projectedDistance;
            minDistFaceLabel = containingCell[index];
        }
    }

    return minDistFaceLabel;
}

label frontSmoother::containingEdge(const label& pointLabel, const label& faceLabel, const triSurfaceFront& front, const fvMesh& mesh) const
{
    label edgeWithFrontVertex = -1;
     
    // Check if the point coincides with an edge
    const auto& frontVertex = front.localPoints()[pointLabel];
    const auto& faceEdges = mesh.faceEdges(faceLabel);
    const auto& meshEdges = mesh.edges();
    const auto& meshVertices = mesh.points();
    
    vector edgeDirection{0,0,0};
    vector e0FrontVertexDirection{0,0,0};

    forAll(faceEdges, index)
    {
        const auto& faceEdge = meshEdges[faceEdges[index]];

        edgeDirection = meshVertices[faceEdge[1]] - meshVertices[faceEdge[0]];
        edgeDirection /= mag(edgeDirection);
        e0FrontVertexDirection = frontVertex - meshVertices[faceEdge[0]];
        e0FrontVertexDirection /= (mag(e0FrontVertexDirection) + SMALL);

        if (mag(1.0 - (edgeDirection&e0FrontVertexDirection)) < SMALL)
        {
            edgeWithFrontVertex = faceEdges[index];
            break;
        }
    }

    return edgeWithFrontVertex;
}

bool frontSmoother::boundaryFacesAreCoplanar(const label& edgeLabel, const triSurfaceFront& front, const fvMesh& mesh) const

{
    const auto& edgeFaces = mesh.edgeFaces(edgeLabel);
    auto nInternalFaces = mesh.nInternalFaces();
    vectorTuple boundaryFaceNormals{vector{0,0,0}, vector{0,0,0}};
    label counter = 0;

    // Assumption: each boundary edge of the fvMesh is shared by
    // exactly two boundary faces.
    // They are coplanar if their normals are parallel
    forAll(edgeFaces, index)
    {
        const auto& faceLabel = edgeFaces[index];

        if (faceLabel < nInternalFaces)
        {
            continue;
        }
        else
        {
            boundaryFaceNormals[counter] = mesh.faceAreas()[faceLabel];
            boundaryFaceNormals[counter] /= mag(boundaryFaceNormals[counter]);

            ++counter;
        }
    }

    if (mag(1.0 - (boundaryFaceNormals[0] & boundaryFaceNormals[1])) < SMALL)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
frontSmoother::frontSmoother(const dictionary& configDict)
:
    underrelaxationFactor_{readScalar(configDict.lookup("relaxFactor"))},
    nSweeps_{readLabel(configDict.lookup("nSweeps"))},
    smoothingType_{configDict.lookup("smooth")}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void frontSmoother::smoothEdges(triSurfaceFront& front, const fvMesh& mesh) const
{
    pointField& points = const_cast<pointField&>(front.localPoints());

    auto singlePointBoundaryEdgeLabels = singlePointBoundaryEdges(front);
    auto internalEdgeLabels = internalEdges(front.nInternalEdges(), singlePointBoundaryEdgeLabels);

    for(label sweep=0; sweep < nSweeps_; ++sweep)
    {
        const auto& edges = front.edges();
        
        // Iterate internal Edges
        for(const auto& edgeLabel : internalEdgeLabels)
        {
            const auto& relaxEdge = edges[edgeLabel];

            vectorTuple Ai{vector{0,0,0}, vector{0,0,0}};
            Ai[0] = computeA(relaxEdge[0], front);
            Ai[1] = computeA(relaxEdge[1], front);

            auto v = computeV(edgeLabel, front);
            
            auto smoothedPoints = smoothEdge(relaxEdge, front);
            auto dxs0 = underrelaxationFactor_*(smoothedPoints[0] - points[relaxEdge[0]]);
            auto dxs1 = underrelaxationFactor_*(smoothedPoints[1] - points[relaxEdge[1]]);

            auto Afinal = Ai[0] + Ai[1] + (v ^ (dxs0 - dxs1));

            // The paper only specifies this threshold as "a tiny number" (TT)
            if (mag(Afinal) > 1.0e-10)
            {
                auto normal = Afinal / mag(Afinal);

                auto h = -((dxs0 & Ai[0]) + (dxs1 & Ai[1])
                            + (dxs1 & (v ^ dxs0))) / mag(Afinal);

                points[relaxEdge[0]] += dxs0 + h*normal;
                points[relaxEdge[1]] += dxs1 + h*normal;
            }
        }

        // Iterate single point boundary edges
        for(const auto& edgeLabel : singlePointBoundaryEdgeLabels)
        {
            const auto& relaxEdge = edges[edgeLabel];

            vectorTuple Ai{vector{0,0,0}, vector{0,0,0}};
            Ai[0] = computeA(relaxEdge[0], front);
            Ai[1] = computeA(relaxEdge[1], front);

            auto v = computeV(edgeLabel, front);
            
            auto smoothedPoints = smoothEdge(relaxEdge, front);

            vectorTuple dxs{vector{0,0,0}, vector{0,0,0}};
            dxs[0] = underrelaxationFactor_*(smoothedPoints[0] - points[relaxEdge[0]]);
            dxs[1] = underrelaxationFactor_*(smoothedPoints[1] - points[relaxEdge[1]]);
            
            // Determine which point of relaxEdge is on boundary
            auto boundaryPointIndex = determineBoundaryPoint(relaxEdge, front);

            // Correct displacement of boundary point so it is aligned with
            // face/edge
            auto moveRestriction = movementRestriction(relaxEdge[boundaryPointIndex], front, mesh);
            dxs[boundaryPointIndex] = moveRestriction & dxs[boundaryPointIndex];

            auto Afinal = Ai[0] + Ai[1] + (v ^ (dxs[0] - dxs[1]));

            // The paper only specifies this threshold as "a tiny number" (TT)
            if (mag(Afinal) > 1.0e-10)
            {
                vectorTuple normals{vector{Afinal / mag(Afinal)}, vector{Afinal / mag(Afinal)}};
                normals[boundaryPointIndex] = moveRestriction & normals[boundaryPointIndex];

                // In case both normals are not parallel, solution of a
                // quadratic equation to obtain h is required
                scalar h = 0;
                auto p1 = normals[1] & (v ^ normals[0]);
                auto p2 = (normals[0]&Ai[0]) + (normals[1]&Ai[1])
                          + (normals[1]&(v^dxs[0])) + (dxs[1]&(v^normals[0]));
                auto p3 = (dxs[0]&Ai[0]) + (dxs[1]&Ai[1]) + (dxs[1]&(v^dxs[0]));

                if (p1 > SMALL)
                {
                    // TODO: determine which sign makes sense for the sqrt term
                    // only '+' works, figure out why
                    h = (-p2 + Foam::sqrt(p2*p2 - 4*p1*p3)) / (2*p1);
                }
                else
                {
                    h = -p3/p2;
                }

                points[relaxEdge[0]] += dxs[0] + h*normals[0];
                points[relaxEdge[1]] += dxs[1] + h*normals[1];
            }
        }

        // FIXME: for now do not smooth boundary edges with both points on
        // the boundary. As so far from the simulations it seems to be
        // sufficient to smooth the one point boundary edges
    }

    updateGlobalPoints(front);
}

void frontSmoother::smoothPoints(triSurfaceFront& front, const fvMesh& mesh) const
{
    pointField& points = const_cast<pointField&>(front.localPoints());

    for (label sweep = 0; sweep < nSweeps_; ++sweep)
    {
        const auto& boundaryPoints = front.boundaryPoints();
        label boundaryPointIndex = 0;

        forAll(points, pointLabel)
        {
            auto sumAi = computeA(pointLabel, front);
            // Compute normal at point as area averaged normal of all connected
            // triangles
            auto normal = sumAi / (mag(sumAi) + SMALL);
            auto barycentre = starBarycentre(pointLabel, front);
            auto dxs = barycentre - points[pointLabel];

            if (boundaryPoints.size() > 0 && pointLabel == boundaryPoints[boundaryPointIndex])
            {
                auto moveRestriction = movementRestriction(pointLabel, front, mesh);

                // Project to boundary plane/edge
                normal = moveRestriction & normal;
                normal /= (mag(normal) + SMALL);
                dxs = moveRestriction & dxs;

                points[pointLabel] += dxs - (dxs&sumAi)/((normal&sumAi) + SMALL)*normal;
                ++boundaryPointIndex;
            }
            else
            {
                points[pointLabel] += dxs - (dxs & normal) * normal;
            }
        }
    }

    updateGlobalPoints(front);
}

void frontSmoother::smoothFront(triSurfaceFront& front, const fvMesh& mesh) const
{
    // TODO: rather implement this as a runtime selectable class hierarchy? (TT)
    if (smoothingType_ == "edges")
    {
        Info << "\nSmoothing edges...\n";
        smoothEdges(front, mesh);
    }
    else if (smoothingType_ == "points")
    {
        Info << "\nSmoothing points...\n";
        smoothPoints(front, mesh);
    }
    else if (smoothingType_ == "pointsAndEdges")
    {
        Info << "\nSmoothing points and edges...\n";
        smoothPoints(front, mesh);
        smoothEdges(front, mesh);
    }
    else
    {
        FatalErrorIn (
            "frontSmoother::smoothFront(front, mesh)"
        )   << "Unknown algorithm type "
            << smoothingType_
            << " for the keyword 'smooth'.\n"
            << "Valid options are 'points', 'edges' or 'pointsAndEdges'."
            << exit(FatalError);
    }

    // IMPORTANT: front smoothing invalidates demand driven front geometry data
    // e.g. face normals
    // --> delete them
    front.clearGeom();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
