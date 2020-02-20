/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include <cmath>

#include "lentFrontSmoothingTest.H"

#include "fvCFD.H"
#include "volMesh.H"

#include "errorMetrics.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
bool lentFrontSmoothingTest::normalsAreConsistent()
{
    return normalConsistencyPtr_->makeFrontNormalsConsistent(frontRef());
}


bool lentFrontSmoothingTest::boundaryPointsRemainOnBoundary() const
{
    // NOTE: for simplicity this check assumes that the parts of  boundary
    // in contact
    // with the front are planar and that their normals are parallel to one
    // of the coordinate axes
    bool bPointsRemainedOnBoundary = true;

    const auto& boundingBox = mesh().bounds();
    auto projector = twoDProjector();
    auto emptyDirection = (Identity<scalar>{} - projector)&vector{1, 1, 1};

    const auto& front = frontRef();
    const auto& vertices = front.localPoints();
    const auto& boundaryPoints = front.boundaryPoints();

    const scalar tolerance = 10*SMALL;

    for (const auto& bpLabel : boundaryPoints)
    {
        if
        (
            mag((vertices[bpLabel] - boundingBox.min()) & emptyDirection) > tolerance 
         && mag((vertices[bpLabel] - boundingBox.max()) & emptyDirection) > tolerance
        )
        {
            bPointsRemainedOnBoundary = false;
        }
    }

    return bPointsRemainedOnBoundary;
}


bool lentFrontSmoothingTest::is2DCase() const
{
    return mesh().nGeometricD() < 3;
}


point lentFrontSmoothingTest::edgeLoopCentre(const labelList& loop, const triSurfaceFront& front) const
{
    point centre{0, 0, 0};

    const auto& v = front.localPoints();
    const auto& edges = front.edges();

    for (const auto& edgeLabel : loop)
    {
        // Only use first point of edge to avoid adding the same point
        // twice during iteration (TT)
        centre += v[edges[edgeLabel][0]];
    }

    return centre/loop.size();
}


scalar lentFrontSmoothingTest::tetVolume(const point& base, const point& v0, const point& v1, const point& v2) const
{
    return mag(((v0 - base) ^ (v1 - base)) & (v2 - base))/6.0;
}


scalar lentFrontSmoothingTest::frontVolume() const
{
    scalar volume = 0.0;

    const auto& front = frontRef();
    const auto& v = front.localPoints();
    const auto& triangles = front.localFaces();

    point geometricCentre{0, 0, 0};

    for (const auto& vertex : v)
    {
        geometricCentre += vertex;
    }
    geometricCentre /= v.size();

    for (const auto& tria : triangles)
    {
        volume += tetVolume(v[tria[0]], v[tria[1]], v[tria[2]], geometricCentre);
    }

    // Add the missing volume of the open front in a 2D case
    if (is2DCase())
    {
        const auto& edgeLoops = front.edgeLoops();
        const auto& edges = front.edges();

        for (const auto& eLoop : edgeLoops)
        {
            auto loopCentre = edgeLoopCentre(eLoop, front);

            for (const auto& eLabel : eLoop)
            {
                const auto& anEdge = edges[eLabel];

                volume += tetVolume(geometricCentre, v[anEdge[0]], v[anEdge[1]], loopCentre);
            }
        }
    } 

    return volume;
}


point lentFrontSmoothingTest::computeFrontCentreOfGravity() const
{
    point cog{0, 0, 0,};

    const auto& front = frontRef();
    const auto& triangles = front.localFaces();
    const auto& v = front.localPoints();

    point geometricCentre{0, 0, 0};

    for (const auto& vertex : v)
    {
        geometricCentre += vertex;
    }
    geometricCentre /= v.size();

    // Inner part: yielding centre of gravity for a closed surface
    for (const auto& tria : triangles)
    {
        cog += 0.25*tetVolume(v[tria[0]], v[tria[1]], v[tria[2]], geometricCentre)
                *(v[tria[0]] + v[tria[1]], v[tria[2]], geometricCentre);
    }
    
    // Boundary part: missing contribution in case of open surfaces
    if (is2DCase())
    {
        const auto& edgeLoops = front.edgeLoops();
        const auto& edges = front.edges();

        for (const auto& eLoop : edgeLoops)
        {
            auto loopCentre = edgeLoopCentre(eLoop, front);

            for (const auto& eLabel : eLoop)
            {
                const auto& anEdge = edges[eLabel];

                cog += 0.25*tetVolume(geometricCentre, v[anEdge[0]], v[anEdge[1]], loopCentre)
                    *(geometricCentre + v[anEdge[0]] + v[anEdge[1]] + loopCentre);
            }

        }
    }

    return cog/frontVolume();
}


tensor lentFrontSmoothingTest::twoDProjector() const
{
    // geometricD() returns a vector in which eahc component assumes
    // either "1" (valid) oder "-1" (invalid)
    auto inValidDimensions = -0.5*(vector{mesh().geometricD()} - vector{1, 1, 1});
    auto projector = Identity<scalar>{} - inValidDimensions*inValidDimensions;

    return projector;
}


void lentFrontSmoothingTest::addFrontNormalNoise()
{
    auto& front = frontRef();
    auto& vertices = const_cast<pointField&>(front.points());
    const auto& pointNormals = front.pointNormals();

    // For 2D cases it must be ensured that the boundary points are
    // not moved in the empty direction. Thus use projection to remove
    // this component from the noise.
    auto projector = twoDProjector();

    forAll(vertices, I)
    {
        auto normalProjector = pointNormals[I]*pointNormals[I];
        vertices[I] += projector&(normalProjector&noiseGen_.noise<vector>(magFrontNoise_));
    }
}


void lentFrontSmoothingTest::computeSurfaceDeviation()
{
    auto& surfaceDeviation = surfaceDeviationTmp_.ref();
    const auto& surface = surfaceRef();
    const auto& vertices = frontRef().localPoints();
    
    // Use the cubic root of the volume as characteristic length
    // for normalization
    auto characteristicLength = std::cbrt(originalVolume_);

    forAll(vertices, I)
    {
        surfaceDeviation[I] = surface.signedDistance(vertices[I])/characteristicLength;
    }
}


void lentFrontSmoothingTest::randomSetup()
{
    setupFrontFromSurface(correctFront_);

    normalConsistencyPtr_->makeFrontNormalsConsistent(frontRef());

    // Fields associated with the front are not automatically resized,
    // thus do it manually here after the front setup
    surfaceDeviationTmp_.ref().resize(frontRef().localPoints().size());
}


void lentFrontSmoothingTest::perturbInputFields()
{
    if (correctFront_ && magFrontNoise_ > SMALL)
    {
        const auto& surface = surfaceRef();
        auto& front = frontRef();

        surface.moveFrontToSurface(front);
        
        addFrontNormalNoise();
       
        // Clear buffered point dependent data
        front.clearGeom();
    }

    originalVolume_ = frontVolume();
    originalCentreOfGravity_ = computeFrontCentreOfGravity();
}


void lentFrontSmoothingTest::computeApproximatedFields()
{
    frontSmootherPtr_->smoothFront(frontRef(), mesh());
}

void lentFrontSmoothingTest::evaluateMetrics()
{
    auto newVolume = frontVolume();
    auto newCentre = computeFrontCentreOfGravity();
    scalar boundaryPointsOK = 0.0;
    scalar consistentNormals = 0.0;

    if (boundaryPointsRemainOnBoundary())
    {
        boundaryPointsOK = 1.0;
    }

    if (normalsAreConsistent())
    {
        consistentNormals = 1.0;
    }

    addMeasure("relative_volume_error", mag(newVolume - originalVolume_)/originalVolume_);
    // NOTE: normalize the centre shift with a characteristic length:
    // here: cubic root of volume
    addMeasure("mag_centre_shift", mag(originalCentreOfGravity_ - newCentre)/std::cbrt(originalVolume_));
    addMeasure("centre_shift", (newCentre - originalCentreOfGravity_));
    addMeasure("boundary_points_ok", boundaryPointsOK);
    addMeasure("consistent_normals", consistentNormals);

    computeSurfaceDeviation();
    errorMetrics surfaceDeviationMetrics{surfaceDeviationTmp_.ref()};

    addMeasure("L_1_position", surfaceDeviationMetrics.arithmeticMeanError());
    addMeasure("L_2_position", surfaceDeviationMetrics.quadraticMeanError());
    addMeasure("L_inf_position", surfaceDeviationMetrics.maximumError());

    frontRef().write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentFrontSmoothingTest::lentFrontSmoothingTest(const fvMesh& mesh, triSurfaceFront& front)
:
    lentSubalgorithmTest{mesh, front},
    frontSmootherPtr_{},
    normalConsistencyPtr_{},
    surfaceDeviationTmp_{}
{
    frontSmootherPtr_ = std::unique_ptr<frontSmoother>
                        {
                            new frontSmoother{lentDict().subDict("frontSmootherTest")}
                        };

    normalConsistencyPtr_ = std::unique_ptr<analyticalSurfaceNormalConsistency>
                                {
                                    new analyticalSurfaceNormalConsistency{surfaceRef()}
                                };

    correctFront_ = testDict().get<Switch>("correctFront");
    magFrontNoise_ = testDict().get<scalar>("magFrontNoise");

    surfaceDeviationTmp_ = tmp<triSurfaceFrontPointScalarField>{
                new triSurfaceFrontPointScalarField
                {
                    IOobject{
                        "surface_deviation", 
                        mesh.time().timeName(), 
                        front,
                        IOobject::NO_READ, 
                        IOobject::AUTO_WRITE
                    }, 
                    frontRef(), 
                    dimensionedScalar{
                        "zero", 
                        dimless, 
                        0.0
                    }
                }
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
