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
    Foam::parabolaFitting2DNormalCalculator

SourceFiles
    parabolaFitting2DNormalCalculator.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Compute the front vertex normals from a parabola fitted at each vertex.
    The parameters of the parabola f(x,y)=ax^2 + by^2 + cxy + dx + ey
    are found by solving a least-squares problem considering all vertices
    of the two-ring neighbourhood.

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

#include "parabolaFitting2DNormalCalculator.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(parabolaFitting2DNormalCalculator, 0);
    addToRunTimeSelectionTable(frontVertexNormalCalculator, parabolaFitting2DNormalCalculator, Dictionary);

// * * * * * * * * * * * *  Private member functions * * * * * * * * * * * * //
void parabolaFitting2DNormalCalculator::initializeVertexNormal(const fvMesh& mesh, const triSurfaceFront& front) const
{
    const Time& runTime = mesh.time();  

    vertexNormalTmp_ = tmp<triSurfaceFrontPointVectorField>
    (
        new triSurfaceFrontPointVectorField
        (
            IOobject(
                "frontNormals", 
                runTime.timeName(), 
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
        )
    );
}

void parabolaFitting2DNormalCalculator::computeProjectionData(const fvMesh& mesh) const
{
    if (projector_ == tensor(Identity<scalar>{}))
    {
        const auto& patches = mesh.boundary();

        forAll(patches, I)
        {
            if (patches[I].type() == "empty")
            {
                // FIXME: works NOT as intended. Find alternative to detect
                // the empty direction
                //
                //auto emptyDirection = patches[I].Sf()[0];
                //emptyDirection /= mag(emptyDirection);
                vector emptyDirection{0, 0, 1};

                //----------------------------------------
                emptyDirection_ = emptyDirection;

                // Distinguish between the three possibilities for the empty direction
                if(mag(emptyDirection_[0]) > SMALL)
                {
                    xRef_ = 1;
                    yRef_ = 2;
                    zRef_ = 0;
                    heightDirection_ = vector{0, 0, 1};
                }
                else if (mag(emptyDirection_[1]) > SMALL)
                {
                    xRef_ = 0;
                    yRef_ = 2;
                    zRef_ = 1;
                    heightDirection_ = vector{0, 0, 1};
                }
                else if (mag(emptyDirection_[2]) > SMALL)
                {
                    xRef_ = 0;
                    yRef_ = 1;
                    zRef_ = 2;
                    heightDirection_ = vector{0, 1, 0};
                }

                break;
            }
        }
    }
}

quaternion parabolaFitting2DNormalCalculator::rotationQuaternion(const triSurfaceFront& front, const label& pointLabel) const
{
    auto approxNormal = front.pointNormals()[pointLabel];

    // transform to reference system
    approxNormal = projector_&approxNormal;
    approxNormal /= mag(approxNormal) + SMALL;

    auto rotationAxis = approxNormal^heightDirection_;

    if (mag(rotationAxis) > SMALL)
    {
        rotationAxis /= mag(rotationAxis);
    }
    else
    {
        rotationAxis = emptyDirection_;
    }

    return quaternion{rotationAxis, heightDirection_&approxNormal, true};
}

vectorXd parabolaFitting2DNormalCalculator::fitParabola(const std::vector<Foam::vector>& points) const
{
    // This function fits the parabola
    // z(x,y) = a0x^2 + a1y^2 + a2xy + a3x + a4y
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
    vertexNormalTmp_{},
    xRef_{0.0},
    yRef_{0.0},
    zRef_{0.0},
    emptyDirection_{0,0,0},
    heightDirection_{0,0,0}
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
tmp<triSurfaceFrontPointVectorField> parabolaFitting2DNormalCalculator::vertexNormals(const fvMesh& mesh, const triSurfaceFront& front) const
{
    computeProjectionData(mesh);

    if (vertexNormalTmp_.empty())
    {
        initializeVertexNormal(mesh, front);
    }

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
