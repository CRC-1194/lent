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
    Foam::parabolaFitting2DNormalCalculator

SourceFiles
    parabolaFitting2DNormalCalculator.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)

Description

    Compute the front vertex normals from a parabola fitted at each vertex.
    The parameters of the parabola f(x,y)=ax^2 + bx
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

#include "parabolaFitting2DNormalCalculator.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(parabolaFitting2DNormalCalculator, 0);
    addToRunTimeSelectionTable(frontVertexNormalCalculator, parabolaFitting2DNormalCalculator, Dictionary);

// * * * * * * * * * * * *  Private member functions * * * * * * * * * * * * //
void parabolaFitting2DNormalCalculator::computeProjectionData(const fvMesh& mesh) const
{
    if (projector_ == tensor(Identity<scalar>{}))
    {
        // Derive empty direction from the mesh bounding
        // box. This assumes that mesh has the least expanse in the
        // empty direction
        const auto& bb = mesh.bounds().span();

        if (bb[0] < bb[1] && bb[0] < bb[2])
        {
            // x is empty direction
            emptyDirection_ = vector{1, 0, 0};
            heightDirection_ = vector{0, 0, 1};
            xRef_ = 1;
            yRef_ = 2;
            zRef_ = 0;
        }
        else if (bb[1] < bb[0] && bb[1] < bb[2])
        {
            // y is empty direction
            emptyDirection_ = vector{0, 1, 0};
            heightDirection_ = vector{0, 0, 1};
            xRef_ = 0;
            yRef_ = 2;
            zRef_ = 1;
        }
        else
        {
            // z is empty direction
            emptyDirection_ = vector{0, 0, 1};
            heightDirection_ = vector{0, 1, 0};
            xRef_ = 0;
            yRef_ = 1;
            zRef_ = 2;
        }

        projector_ = Identity<scalar>{} - emptyDirection_*emptyDirection_;
    }
}

quaternion parabolaFitting2DNormalCalculator::rotationQuaternion(const triSurfaceFront& front, const label& pointLabel) const
{
    auto approxNormal = front.pointNormals()[pointLabel];

    // transform to reference system
    approxNormal[zRef_] = 0.0;
    approxNormal /= (mag(approxNormal) + SMALL);

    vector rotationAxis{0,0,0};

    // Rotation axis and the corresponding angle have to be determined
    // carefully in order to not compromise the accuracy of the solution.
    if (approxNormal[xRef_] > 0.0)
    {
        rotationAxis = emptyDirection_;
    }
    else
    {
        rotationAxis = -1.0*emptyDirection_;
    }

    return quaternion{rotationAxis, heightDirection_&approxNormal, true};
}

vectorXd parabolaFitting2DNormalCalculator::fitParabola(const std::vector<Foam::vector>& points) const
{
    // This function fits the parabola
    // y(x) = ax^2 + bx
    // to the given point set.The aboslute term is missing so the parabola
    // always includes the point we are actually fitting to
    matrixXd A{points.size(), 2};
    vectorXd b{points.size()};

    label index = 0;
    for (const auto& aPoint : points)
    {
        A(index, 0) = aPoint[xRef_]*aPoint[xRef_];
        A(index, 1) = aPoint[xRef_];

        b(index) = aPoint[yRef_];

        ++index;
    }

    return A.colPivHouseholderQr().solve(b);
}

vector parabolaFitting2DNormalCalculator::computeNormal(const vectorXd& coefficients, const quaternion& rotation) const
{
    vector localNormal{0, 0, 0};
    localNormal[xRef_] = -1.0*coefficients(1);
    localNormal[yRef_] = 1.0;

    localNormal /= mag(localNormal) + SMALL;

    return rotation.invTransform(localNormal);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
parabolaFitting2DNormalCalculator::parabolaFitting2DNormalCalculator(const dictionary& configDict)
:
    frontVertexNormalCalculator{configDict},
    projector_{Identity<scalar>{}},
    xRef_{0},
    yRef_{0},
    zRef_{0},
    emptyDirection_{0,0,0},
    heightDirection_{0,0,0}
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<triSurfaceFrontPointVectorField> parabolaFitting2DNormalCalculator::vertexNormals(const fvMesh& mesh, const triSurfaceFront& front) const
{
    computeProjectionData(mesh);

    tmp<triSurfaceFrontPointVectorField> vertexNormalTmp_
    {
        new triSurfaceFrontPointVectorField
        {
            IOobject(
                "frontNormals", 
                mesh.time().timeName(), 
                front,
                IOobject::NO_READ, 
                IOobject::AUTO_WRITE
            ), 
            front, 
            dimensionedVector(
                "zero", 
                dimless, 
                vector(0.0,0.0,0.0)
            )
        }
    };

    auto& normals = vertexNormalTmp_.ref();

    const auto& twoRingNeighbours = front.twoRingNeighbours();
    const auto& points = front.localPoints();

    forAll(normals, I)
    {
        const auto& pointSet = twoRingNeighbours[I];

        auto rotationQ = rotationQuaternion(front, I);

        std::vector<Foam::vector> localCoordinatePoints{};

        for (const auto& pointLabel : pointSet)
        {
            auto localVector = projector_&(points[pointLabel] - points[I]);
            localVector = rotationQ.transform(localVector);
            localCoordinatePoints.push_back(localVector);
        }

        auto parabolaCoefficients = fitParabola(localCoordinatePoints);

        auto normal = computeNormal(parabolaCoefficients, rotationQ); 

        // Finally, check orientation
        if ((normal & front.pointNormals()[I]) < 0.0)
        {
            normal *= -1.0;
        }

        normals[I] = normal/mag(normal);
    }

    return vertexNormalTmp_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
