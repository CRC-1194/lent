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

#include "errorMetrics.H"
#include "lentVertexNormalCalculatorTest.H"

#include <cmath>

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
tensor lentVertexNormalCalculatorTest::twoDProjector() const
{
    // geometricD() returns a vector in which eahc component assumes
    // either "1" (valid) oder "-1" (invalid)
    auto inValidDimensions = -0.5*(vector{mesh().geometricD()} - vector{1, 1, 1});
    auto projector = Identity<scalar>{} - inValidDimensions*inValidDimensions;

    return projector;
}

void lentVertexNormalCalculatorTest::addFrontNormalNoise()
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void lentVertexNormalCalculatorTest::randomSetup()
{
    setupFrontFromSurface(correctFront_);

    // Ensure consistency of the front triangle normals, otherwise the
    // check for inconsistent normals is senseless
    normalConsistencyPtr_->makeFrontNormalsConsistent(frontRef());
}

void lentVertexNormalCalculatorTest::perturbInputFields()
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
}

void lentVertexNormalCalculatorTest::computeApproximatedFields()
{
    const auto& normalCalculator = normalCalculatorTmp_.ref();

    approximateNormalsTmp_ = normalCalculator.vertexNormals(mesh(), frontRef());
}

void lentVertexNormalCalculatorTest::evaluateMetrics()
{
    const auto& approximateNormals = approximateNormalsTmp_.ref();
    const auto& surface = surfaceRef();
    const auto& vertices = frontRef().localPoints();

    normalDeviationField_.resize(approximateNormals.size());

    vector exactNormal{0,0,0};

    bool normalsHaveUnitLength = true;
    bool normalsAreConsistent = true;

    forAll(vertices, I)
    {
        if (mag(mag(approximateNormals[I]) - 1.0) > 10*SMALL)
        {
            normalsHaveUnitLength = false;
        }

        exactNormal = surface.normalToPoint(vertices[I]);

        if ((exactNormal & approximateNormals[I]) < 0.0)
        {
            normalsAreConsistent = false;
        }

        auto dotProduct = exactNormal & approximateNormals[I];

        // Cut off over-/undershoots of the dot product to avoid exceptions
        // in std::acos. 
        if (mag(dotProduct) > 1.0)
        {
            dotProduct = 1.0*sign(dotProduct);
        }

        normalDeviationField_[I] = std::acos(dotProduct);
    }

    errorMetrics normalDeviationMetrics{normalDeviationField_};

    addMeasure("unit_length", scalar(normalsHaveUnitLength));
    addMeasure("consistent_normals", scalar(normalsAreConsistent));

    addMeasure("L_1_normal_error", normalDeviationMetrics.arithmeticMeanError());
    addMeasure("L_2_normal_error", normalDeviationMetrics.quadraticMeanError());
    addMeasure("L_inf_normal_error", normalDeviationMetrics.maximumError());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentVertexNormalCalculatorTest::lentVertexNormalCalculatorTest(const fvMesh& mesh, triSurfaceFront& front)
:
    lentSubalgorithmTest{mesh, front},
    normalCalculatorTmp_{},
    normalConsistencyPtr_{},
    approximateNormalsTmp_{},
    normalDeviationField_{}
{
    normalCalculatorTmp_ = tmp<frontVertexNormalCalculator>{
                                frontVertexNormalCalculator::New
                                (
                                    lentDict().subDict("normalCalculator")
                                )
                            };

    normalConsistencyPtr_ = std::unique_ptr<analyticalSurfaceNormalConsistency>
                                {
                                    new analyticalSurfaceNormalConsistency{surfaceRef()}
                                };

    correctFront_ = Switch{testDict().lookup("correctFront")};
    magFrontNoise_ = readScalar(testDict().lookup("magFrontNoise"));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
