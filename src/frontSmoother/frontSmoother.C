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

#include "frontSmoother.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void frontSmoother::updateLocalPoints(triSurfaceFront& front) const
{
    pointField& localPoints = const_cast<pointField&>(front.localPoints());
    const auto& frontPoints = front.points();
    const labelList& globalToLocal = front.meshPoints();

    forAll(globalToLocal, pointI)
    {
        localPoints[pointI] = frontPoints[globalToLocal[pointI]];
    }
}

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

FixedList<vector,2> frontSmoother::smoothEdge(const edge& relaxEdge, const triSurfaceFront& front) const
{
    FixedList<vector,2> smoothedPoints{vector{0,0,0}, vector{0,0,0}};

    smoothedPoints[0] = starBarycentre(relaxEdge[0], front);
    smoothedPoints[1] = starBarycentre(relaxEdge[1], front);

    return smoothedPoints;
}

vector frontSmoother::boundaryNormal(const label& pointLabel, const triSurfaceFront& front, const fvMesh& mesh) const
{
    vector normal{0,0,0};

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
    label minDistFaceLabel = 0;

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

    // TODO: ensure these normals are oriented outwards of the mesh
    // --> think again, for projection the sign of the normal is irrelevant
    normal = mesh.faceAreas()[minDistFaceLabel];

    return normal/mag(normal);
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


void frontSmoother::ensureNormalConsistency(triSurfaceFront& front, const edge& relaxEdge, const vector& normal) const
{
    const auto& pointToFace = front.pointFaces();
    const auto& points = front.localPoints();
    const auto& faces = front.localFaces();

    forAll(relaxEdge, index)
    {
        const auto& pointLabel = relaxEdge[index];
        const auto& connectedFaces = pointToFace[pointLabel];

        forAll(connectedFaces, I)
        {
            const auto& aFace = faces[connectedFaces[I]];

            if ((normal & aFace.normal(points)) < 0)
            {
                labelledTri& modFace = const_cast<labelledTri&>(aFace);
                modFace.flip();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
frontSmoother::frontSmoother(const dictionary& configDict)
:
    underrelaxationFactor_{readScalar(configDict.lookup("relaxFactor"))},
    nSweeps_{readLabel(configDict.lookup("nSweeps"))}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void frontSmoother::smoothEdges(triSurfaceFront& front) const
{
    updateLocalPoints(front);
    pointField& points = const_cast<pointField&>(front.localPoints());

    // Extension for 2D: edges in list returned by edges() are sorted such that
    // all internal edges come before the boundary ones. This means edges with
    // a label X >= nInternalEdges() are out of the game, removing edges where
    // both points are boundary points. However, we still have to eliminate
    // edges where one point is a boundary point.
    // Best idea so far: set a flag if case is 2D and use the point coordinates
    // of the edge to detect whether it is a full interior edge. Other
    // approaches are probably to costly (to many search operations) with the
    // available information of the triSurfaceFront class. 
    // More elegant/robust/general approach yet more involved would be to extend
    // the triSurfaceFront class to obtain and store this information. 

    for(label sweep=0; sweep < nSweeps_; ++sweep)
    {
        const auto& edges = front.edges();

        forAll(edges, edgeLabel)
        {
            const auto& relaxEdge = edges[edgeLabel];

            FixedList<vector,2> Ai{vector{0,0,0}, vector{0,0,0}};
            Ai[0] = computeA(relaxEdge[0], front);
            Ai[1] = computeA(relaxEdge[1], front);

            auto v = computeV(edgeLabel, front);
            
            auto smoothedPoints = smoothEdge(relaxEdge, front);
            auto dxs0 = underrelaxationFactor_*(smoothedPoints[0] - points[relaxEdge[0]]);
            auto dxs1 = underrelaxationFactor_*(smoothedPoints[1] - points[relaxEdge[1]]);

            auto Afinal = Ai[0] + Ai[1] + (v ^ (dxs0 - dxs1));

            // TODO: test different thresholds (TT)
            if (mag(Afinal) > 1.0e-10)
            {
                auto normal = Afinal / mag(Afinal);

                auto h = ((dxs0 & Ai[0]) + (dxs1 & Ai[1])
                            + (dxs1 & (v ^ dxs0))) / mag(Afinal);

                points[relaxEdge[0]] += dxs0 - h*normal;
                points[relaxEdge[1]] += dxs1 - h*normal;
            }

            //ensureNormalConsistency(front, relaxEdge, Ai[0] + Ai[1]);
        }
    }

    updateGlobalPoints(front);
}

void frontSmoother::smoothPoints(triSurfaceFront& front, const fvMesh& mesh) const
{
    updateLocalPoints(front);
    pointField& points = const_cast<pointField&>(front.localPoints());
    
    for (label sweep = 0; sweep < nSweeps_; sweep++)
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

            if (pointLabel == boundaryPoints[boundaryPointIndex])
            {
                auto bNormal = boundaryNormal(pointLabel, front, mesh);

                // Project to boundary plane
                normal = (Identity<scalar>{} - bNormal*bNormal) & normal;
                normal /= (mag(normal) + SMALL);
                dxs = (Identity<scalar>{} - bNormal*bNormal) & dxs;

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
