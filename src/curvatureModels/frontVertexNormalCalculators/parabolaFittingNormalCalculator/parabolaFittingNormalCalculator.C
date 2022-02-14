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
    Foam::parabolaFittingNormalCalculator

SourceFiles
    parabolaFittingNormalCalculator.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)

Description

    Compute the front vertex normals from a parabola fitted at each vertex.
    The parameters of the parabola f(x,y)=ax^2 + by^2 + cxy + dx + ey
    are found by solving a least-squares problem considering all vertices
    of the two-ring neighbourhood.

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

#include "parabolaFittingNormalCalculator.H"
#include "addToRunTimeSelectionTable.H"

#include <chrono>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(parabolaFittingNormalCalculator, 0);
    addToRunTimeSelectionTable(frontVertexNormalCalculator, parabolaFittingNormalCalculator, Dictionary);

// * * * * * * * * * * * *  Private member functions * * * * * * * * * * * * //
quaternion parabolaFittingNormalCalculator::rotationQuaternion(const triSurfaceFront& front, const label& pointLabel) const
{
    const auto& approxNormal = front.pointNormals()[pointLabel];
    const vector zAxis{0.0, 0.0, 1.0};

    auto rotationAxis = (zAxis ^ approxNormal);

    // Avoid null vector for rotation
    if (mag(rotationAxis) < 1.0e-13)
    {
        rotationAxis = vector{0.0, -1.0, 0.0};
    }
    else
    {
        rotationAxis /= (-1.0*mag(rotationAxis));
    }

    return quaternion{rotationAxis, (zAxis & approxNormal), true};
}

vectorXd parabolaFittingNormalCalculator::fitParabola(const std::vector<Foam::vector>& points) const
{
    // This function fits the parabola
    // z(x,y) = a0x^2 + a1y^2 + a2xy + a3x + a4y
    // to the given point set.The aboslute term is missing so the parabola
    // always includes the point we are actually fitting to
    matrixXd A{points.size(), 5};
    vectorXd b{points.size()};

    label index = 0;
    for (const auto& aPoint : points)
    {
        A(index, 0) = aPoint[0]*aPoint[0];
        A(index, 1) = aPoint[1]*aPoint[1];
        A(index, 2) = aPoint[0]*aPoint[1];
        A(index, 3) = aPoint[0];
        A(index, 4) = aPoint[1];

        b(index) = aPoint[2];

        ++index;
    }

    return A.colPivHouseholderQr().solve(b);
}

vector parabolaFittingNormalCalculator::computeNormal(const vectorXd& coefficients, const quaternion& rotation) const
{
    // The normal is computed as the cross product of the tangent vectors
    // obtained by derivation of the parabola in x- an y-direction.
    // Finally, the normal vector is transformed from the local coordinate
    // system back to the global one
    vector tangentX{1.0, 0.0, coefficients(3)};
    vector tangentY{0.0, 1.0, coefficients(4)};

    auto localNormal = tangentX ^ tangentY;
    localNormal /= (mag(localNormal) + SMALL);

    return rotation.invTransform(localNormal);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
parabolaFittingNormalCalculator::parabolaFittingNormalCalculator(const dictionary& configDict)
:
    frontVertexNormalCalculator{configDict}
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<triSurfaceFrontPointVectorField> parabolaFittingNormalCalculator::vertexNormals(const fvMesh& mesh, const triSurfaceFront& front) const
{
    auto start = std::chrono::system_clock::now();

    const Time& runTime = mesh.time();  

    tmp<triSurfaceFrontPointVectorField> normalsTmp
    (
        new triSurfaceFrontPointVectorField
        (
            IOobject(
                "frontNormals", 
                runTime.timeName(), 
                front,
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ), 
            front, 
            dimensionedVector(
                "zero", 
                dimless, 
                vector(0.0,0.0,0.0)
            )
        )
    );

    auto& normals = normalsTmp.ref();

    const auto& twoRingNeighbours = front.twoRingNeighbours();
    const auto& points = front.localPoints();

    forAll(normals, I)
    {
        const auto& pointSet = twoRingNeighbours[I];

        auto rotationQ = rotationQuaternion(front, I);

        // Rather than to actually use a local coordinate system created by
        // rotation, the vectors are rotated using the negative rotation axis.
        // The result is the same regardign the vector decomposition,
        // but the latter approach is simpler to implement
        std::vector<Foam::vector> localCoordinatePoints{};

        for (const auto& pointLabel : pointSet)
        {
            auto localVector = rotationQ.transform(points[pointLabel] - points[I]);
            localCoordinatePoints.push_back(localVector);
        }

        auto parabolaCoefficients = fitParabola(localCoordinatePoints);

        auto normal = computeNormal(parabolaCoefficients, rotationQ); 

        // Finally, check orientation
        if ((normal & front.pointNormals()[I]) < 0.0)
        {
            normal *= -1.0;
        }
        
        normals[I] = normal;
    }
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;

    Info << "Computation time normals: " << elapsedSeconds.count() << endl;

    return normalsTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
