/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

Application
    lentConstructFront

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "lentMethod.H"
#include "analyticalSurface.H"
#include "analyticalPlane.H"

#include <utility>

using namespace FrontTracking;

bool differentSign(scalar a, scalar b)
{
    if (a*b < 0.0) return true;
    else return false;
}

bool notInList(label ID, labelList& list)
{
    forAll(list, I)
    {
        if (list[I] == ID)
        {
            return false;
        }
    }

    return true;
}

label pointIndex(label ID, labelList& list)
{
    label result = -1;

    forAll(list, I)
    {
        if (list[I] == ID) result = I;
    }

    return result;
}

point geoCentre(const labelList& pointIDs, const pointField& points)
{
    point centre(0.0, 0.0, 0.0);

    forAll(pointIDs, I)
    {
        centre += points[pointIDs[I]];
    }

    centre /= pointIDs.size();

    return centre;
}

void orderPointsAngle(labelList& pointIDs, const pointField& points,
                     const point& centre)
{
    // This is probably the most naive approach, yet it is computationally
    // cheap
    //
    // Remember: angles are not measured in radian, but as the cosine
    // value of the angle [ cos(angle) = <a,b>/(|a||b|) ]
    // Should work as long as polygons are sufficiently convex
    //
    // TODO: potential source of error identified:
    // Angle measurement using the cosine only works reliably if no angles
    // greater than pi are involved. Thus think about a solution to identify
    // angles greater than pi
    // Approach: use a plane and the signed distance of the points to
    // distinguish angles < pi and angles > pi. The plane is spanned by
    // the vector centre --> first point and the normal at the centre
    vector refEdge;
    vector testEdge;
    scalar minAngle;
    scalar currentAngle;
    label minPosition;

    for(label I = 0; I < pointIDs.size()-1; I++)
    {
        refEdge = points[pointIDs[I]] - centre;
        minAngle = -1.1; // Has to be smaller the the minmal cosine value
        minPosition = I;

        for(label K = I+1; K < pointIDs.size(); K++)
        {
            testEdge = points[pointIDs[K]] - centre;
            currentAngle = refEdge & testEdge /(mag(refEdge) * mag(testEdge));

            if (currentAngle > minAngle)
            {
                minAngle = currentAngle;
                minPosition = K;
            }
        }

        std::swap(pointIDs[I+1], pointIDs[minPosition]);
    }
}

bool insidePrism(const point& intersect, const point& centre,
                 const point& refPoint, const vector& normalToAxis)
{
    bool inside = true;

    vector unitNormal = normalToAxis / mag(normalToAxis);

    // Use projected distance along unitNormal
    if ( ((intersect - centre) & unitNormal) >
            ((refPoint - centre) & unitNormal) )
    {
        inside = false;
    }

    return inside;
}

void orderPointsPlaneLock(labelList& pointIDs, const pointField& points,
                          const point& centre, const vector& refNormal)
{
    // More elaborate approach. Maybe it works, maybe not...
    // TODO: add more detailled description
    // FIXME: bugged, results are worse than the naive angle approach
    scalarListList signedDistances(pointIDs.size());
    List<analyticalPlane> planes(pointIDs.size());

    vector centreToP;

    // Construct plane for each point, spanned by Pi, centre and the normal
    // direction of the analytical surface
    forAll(pointIDs, I)
    {
        point p = points[pointIDs[I]];
        centreToP = p - centre;
        planes[I] = analyticalPlane(p, (centreToP ^ refNormal));
    }

    // Compute signed distance matrix
    forAll(pointIDs, I)
    {
        scalarList dist(pointIDs.size());
        forAll(planes, K)
        {
            dist[K] = planes[K].signedDistance(points[pointIDs[I]]);
        }
        signedDistances[I] = dist;
    }

    // Now order points
    for (label I = 0; I < pointIDs.size() - 1; I++)
    {
        label swapPointID = I;

        for (label K = I+1; K < pointIDs.size(); K++)
        {
            bool doSwap = true;

            // Check if and how the point pair I,K intersect the planes.
            // On this basis it is decided if K is a suitable candidate for
            // swapping
            for (label L = 0; L < pointIDs.size(); L++)
            {
                if (differentSign(signedDistances[I][L], signedDistances[K][L]))
                {
                    point intersect =
                        planes[L].intersection(points[pointIDs[I]], 
                                points[pointIDs[K]]);
                    // Check if intersection is inside the prism (meaning located
                    // between the centre and the refPoint of plane L)
                    // or outside. In the latter case the intersection does
                    // not invalidate the suitability for swapping
                    vector refDirection = refNormal ^ planes[L].normal();
                    
                    if (insidePrism(intersect, centre, points[pointIDs[L]],
                                    refDirection))
                    {
                        doSwap = false;
                    }
                }
            }

            // Cancel here, if point K has passed the above test
            if (doSwap)
            {
                swapPointID = K;
                break;
            }
        }

        std::swap(pointIDs[I+1], pointIDs[swapPointID]);
    }
}

// Order intersection points in such a way that the rotational direction
// aligns with the normal of the analytical surface
void orientToNormal(labelList& pointIDs, const pointField& points,
                    const point& centre, const vector& refNormal)
{
    vector triNormal = (points[pointIDs[0]] - centre) ^
                            (points[pointIDs[1]] - centre);
    if ( (refNormal & triNormal) < 0.0 )
    {
        for (label I = 0; I < pointIDs.size() / 2; I++)
        {
            std::swap(pointIDs[I], pointIDs[pointIDs.size()-I-1]);
        }
    }
}

void createTriangles(List<triFace>& trias, const labelList& pointIDs, 
                     const label& centreID)
{
    for (label I = 0; I < pointIDs.size()-1; I ++)
    {
        trias.append(triFace(centreID, pointIDs[I], pointIDs[I+1]));
    }

    // Manual addition of last element since the labelList is not cyclic
    trias.append(triFace(centreID, pointIDs[pointIDs.size()-1], pointIDs[0]));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

using namespace FrontTracking;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    IOdictionary lentControlDict(
        IOobject(
           "lentSolution",
           "system",
           mesh.time(),
           IOobject::MUST_READ_IF_MODIFIED,
           IOobject::AUTO_WRITE
        )
    );

    // Required for vertex located fields aka pointXFields
    pointMesh pmesh(mesh);

    tmp<analyticalSurface> analyticalSurfaceTmp(
        analyticalSurface::New(lentControlDict.subDict("analyticalSurface"))
    );

    // Step 1a: Calculation of the vertex based signed Distances
    pointScalarField pointSignedDistances
    (
        IOobject(
            "pointSignedDistances",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pmesh,
        dimensionedScalar(
            "zero",
            dimLength,
            0.0
        )
    );

    const pointField& vertices = mesh.points();

    forAll(vertices, I)
    {
        pointSignedDistances[I] =
                analyticalSurfaceTmp->signedDistance(vertices[I]);
    }

    // Step 1b: search edges intersected by analyticalSurface
    //          --> different distance signs for points
    const edgeList& edges = mesh.edges();
    labelList intersectedEdges(0, 0);

    forAll(edges, I)
    {
        edge E = edges[I];
        label pointA = E[0];
        label pointB = E[1];

        // Pure intersection
        if (differentSign(pointSignedDistances[pointA],
                    pointSignedDistances[pointB]))
        {
            intersectedEdges.append(I);
        }
        // Start or end of edge is already located on surface
        // TODO: this has to be solved in another way, since in this case
        // the point belongs to multiple edges --> no unique mapping
        else if (pointSignedDistances[pointA] == 0.0 ||
                 pointSignedDistances[pointB] == 0.0)
        {
            intersectedEdges.append(I);
        }
    }


    // Step 2: compute intersections
    pointField intersections(intersectedEdges.size());

    forAll(intersections, I)
    {
        edge E = edges[intersectedEdges[I]];

        const point& pointA = vertices[E[0]];
        const point& pointB = vertices[E[1]];

        intersections[I] = analyticalSurfaceTmp->intersection(
                pointA, pointB);
    }
    

    // Create list of all intersected cells from intersected edges
    const labelListList& edgeToCell = mesh.edgeCells();
    labelList intersectedCells(0, 0);

    forAll(intersectedEdges, I)
    {
        label edgeID = intersectedEdges[I];

        labelList cells = edgeToCell[edgeID];

        forAll(cells, K)
        {
            if (notInList(cells[K], intersectedCells))
            {
                intersectedCells.append(cells[K]);
            }
        }
    }


    // Now it's time to create the per-cell triangulation of the surface
    const labelListList& cellToEdge = mesh.cellEdges();
    List<triFace> frontTriangles(0);
    
    labelList intersectsPerCell(0);
    point tmp(0.0, 0.0, 0.0);

    forAll(intersectedCells, I)
    {
        labelList edges = cellToEdge[intersectedCells[I]];
        intersectsPerCell.clear();

        // Collect all intersections for cell I
        forAll(edges, K)
        {
            // Remember: the postion of an intersected edge in the list
            // intersectedEdges is also the postition of the corresponding
            // intersection point
            label pointID = pointIndex(edges[K], intersectedEdges);

            if (pointID > -1)
            {
                intersectsPerCell.append(pointID);
            }
        }

        // Compute geometrical centre and move to the surface
        tmp = geoCentre(intersectsPerCell, intersections);

        // Avoid overlapping triangles
        //orderPointsAngle(intersectsPerCell, intersections, tmp);
        orderPointsPlaneLock(intersectsPerCell, intersections, tmp,
                             analyticalSurfaceTmp->normalToPoint(tmp));

        // Projection the centre after ordering should improve stability
        // of ordering
        tmp = analyticalSurfaceTmp->normalProjectionToSurface(tmp);
        intersections.append(tmp);

        // Ensure correct orientation
        vector surfaceNormal = analyticalSurfaceTmp->
                                normalToPoint(tmp);
        orientToNormal(intersectsPerCell, intersections, tmp, surfaceNormal);

        // Finally create the list of triangles
        createTriangles(frontTriangles, intersectsPerCell,
                        intersections.size()-1);
    }

    // Create and write the front...
    triSurface impendingDoom(frontTriangles, intersections);
    impendingDoom.write("impendingDoom.stl");

    forAll(frontTriangles, I)
    {
        triFace tmp = frontTriangles[I];

        point p0 = intersections[tmp[0]];
        point p1 = intersections[tmp[1]];
        point p2 = intersections[tmp[2]];

        point centre = (p0 + p1 + p2) / 3.0;
        vector refNormal = analyticalSurfaceTmp->normalToPoint(centre);
        vector triNormal = (p1-p0) ^ (p2-p0);
        scalar cosine = (refNormal & triNormal) /
                        (mag(refNormal) * mag(triNormal));

        if ((cosine) < 0.9)
        {
            Info << "Fucked up triangles: "
                 << cosine << endl;
        }
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
