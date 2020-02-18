/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, 
     \\/     M anipulation  |                    TU Darmstadt
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

Class
    Foam::linearLeastSquaresIsoPointCalculator

Description
    Calculates points on an iso-surface \Gamma := {x, f(x) = s}, where 
    s is the isovalue using an RBF interpolation. 

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    linearLeastSquaresIsoPointCalculatorI.H
    linearLeastSquaresIsoPointCalculator.C

\*--------------------------------------------------------------------------*/

#include "linearLeastSquaresIsoPointCalculator.H"
#include "volMesh.H"
#include <cassert> 
#include <iostream>
#include <limits>

namespace Foam { namespace FrontTracking {
    
// Constructor
linearLeastSquaresIsoPointCalculator::linearLeastSquaresIsoPointCalculator(const bool weighted)
:
    edgePoints_{},
    edgeLabels_{},
    contourPoints_{},
    pointCellLabels_{},
    weighted_{weighted}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool linearLeastSquaresIsoPointCalculator::isoValueInInterval
(
    const scalar isoValue, 
    const scalar phi0, 
    const scalar phi1
) const 
{
    if (((phi0 <= isoValue) && (isoValue <= phi1)) ||
        ((phi1 <= isoValue) && (isoValue <= phi0)))
        return true; 

    return false; 
}

bool linearLeastSquaresIsoPointCalculator::equalUnderTolerance(
    const scalar s1, 
    const scalar s2, 
    const scalar tol
) const
{
    return (mag(s1 - s2) < tol * max(max(mag(s1), mag(s2)), 1.0));
}

void linearLeastSquaresIsoPointCalculator::calcEdgePoints
(
    const pointScalarField& pointPhi, 
    const scalar isoValue 
) 
{
    const auto& mesh = pointPhi.mesh().mesh(); 
    
    const auto& meshEdges = mesh.edges(); 
    const auto& meshPoints = mesh.points(); 

    // Resizing here in case of dynamic adaptive mesh refinement.
    // TODO: Reserve the capacity in the class constructor. TM.
    edgePoints_.resize(0);
    edgeLabels_.resize(static_cast<unsigned int>(mesh.nEdges()));
    std::fill(edgeLabels_.begin(), edgeLabels_.end(), -1);

    forAll(meshEdges, edgeI)
    {
        const auto& edge = meshEdges[edgeI];   
        const auto e0 = edge[0]; 
        const auto e1 = edge[1]; 

        // If the edge is intersected at the root point.
        if (isoValueInInterval(isoValue, pointPhi[e0], pointPhi[e1]))
        {
            // Compute the function derivatives at edge points. 
            vector eVec = meshPoints[e1] - meshPoints[e0]; 

            // Compute the interface root parameter using a secant guess.
            double s = linearRoot(isoValue, pointPhi[e0], pointPhi[e1]); 

            // Append the edge iso-point and mark the edge as intersected. 
            edgePoints_.push_back(meshPoints[e0]  + s * eVec);
            edgeLabels_[static_cast<unsigned int>(edgeI)] = static_cast<label>(edgePoints_.size() - 1); 
        }
    }
}

void linearLeastSquaresIsoPointCalculator::calcContourPoints
(
    const volScalarField& cellPhi, 
    const pointScalarField& pointPhi, 
    const scalar isoValue 
)
{
    calcEdgePoints(pointPhi, isoValue); 
    
    const auto& mesh = cellPhi.mesh();
    const auto& meshCellEdgesLists = mesh.cellEdges(); 

    pointCellLabels_.resize(static_cast<unsigned int>(mesh.nCells()));
    std::fill(pointCellLabels_.begin(), pointCellLabels_.end(), -1);
    contourPoints_.resize(0);

    forAll(meshCellEdgesLists, cellI)
    {
        const auto& cellEdgeLabels = meshCellEdgesLists[cellI]; 
        label nEdgePoints = 0;
        auto startPoint = vector(0,0,0);
        forAll(cellEdgeLabels, edgeI)  
        {
            const label edgeGl = cellEdgeLabels[edgeI]; 
            if (edgeLabels_[edgeGl] > 0)
            {
                startPoint += edgePoints_[edgeLabels_[edgeGl]];  
                ++nEdgePoints; 
            }
        }

        if (nEdgePoints > 0)
        {
            startPoint /= nEdgePoints;
            auto contourPoint = projectToLeastSquaresPlane(startPoint, cellI, cellPhi, pointPhi);
            contourPoints_.push_back(contourPoint); 
            pointCellLabels_[static_cast<unsigned int>(cellI)] = static_cast<label>(contourPoints_.size() - 1);
        }
    }
}

point linearLeastSquaresIsoPointCalculator::projectToLeastSquaresPlane
(
    point aPoint,
    label cellI,
    const volScalarField& cellPhi,
    const pointScalarField& pointPhi
) const
{
    const auto& cellToPoints = cellPhi.mesh().cellPoints()[cellI];
    const auto& meshVertices = cellPhi.mesh().points();

    // Compute least squares plane
    MatrixX4d A{cellToPoints.size() + 1, 4};
    VectorXd phiPrescribed{cellToPoints.size() + 1};

    // Add cell vertex distances
    forAll(cellToPoints, I)
    {
        const auto vertexI = cellToPoints[I];
        auto v = meshVertices[vertexI];
        A(I,0) = 1.0;
        A(I,1) = v.x();
        A(I,2) = v.y();
        A(I,3) = v.z();

        phiPrescribed(I) = pointPhi[vertexI];
    }

    // Add cell centre distance
    auto c = cellPhi.mesh().C()[cellI];
    A(cellToPoints.size(),0) = 1.0;
    A(cellToPoints.size(),1) = c.x();
    A(cellToPoints.size(),2) = c.y();
    A(cellToPoints.size(),3) = c.z();
    phiPrescribed(cellToPoints.size()) = cellPhi[cellI];

    auto weights = computeWeights(aPoint, cellI, cellPhi, pointPhi);

    VectorXd plane = (weights*A).fullPivHouseholderQr().solve(weights*phiPrescribed);
    vector n{plane(1), plane(2), plane(3)};
    scalar d = plane(0);
    auto nMag = mag(n) + EPSILON; 
    n /= nMag;
    d /= nMag;
    auto lambda = -(d + (aPoint&n));

    return aPoint + lambda*n;
}

linearLeastSquaresIsoPointCalculator::DiagonalMatrixXd linearLeastSquaresIsoPointCalculator::computeWeights
(
    const point aPoint,
    const label cellI,
    const volScalarField& cellPhi,
    const pointScalarField& pointPhi
) const
{
    const auto& mesh = cellPhi.mesh();
    const auto& cellToPoints = mesh.cellPoints()[cellI];
    const auto& vertices = mesh.points();

    DiagonalMatrixXd weights(cellToPoints.size() + 1);

    if (!weighted_)
    {
        weights.setIdentity();
        return weights;
    }
        
    //TODO Question: what is a reasonable support length? (TT)
    //TODO Question: appropriate weighting function? (TT)
    //See here: http://www.nealen.net/projects/mls/asapmls.pdf
    
    // Preliminary: compute support width as diagonal of a hex cell for now
    scalar h = 1; //Foam::sqrt(3.0)*Foam::cbrt(cellPhi.mesh().V()[cellI]);
    forAll(cellToPoints, I)
    {
        auto vI = cellToPoints[I];
        auto d = weightInput(vertices[vI], aPoint, mag(pointPhi[vI]));
        weights.diagonal()(I) = weight(d, h);
    }
    auto c = mesh.C()[cellI];
    auto dc = weightInput(c, aPoint, mag(cellPhi[cellI]));
    weights.diagonal()(cellToPoints.size()) = weight(dc, h);

    return weights;
}

scalar linearLeastSquaresIsoPointCalculator::weight(const scalar d, const scalar) const
{
    return 1.0/(d*d + EPSILON);
}

scalar linearLeastSquaresIsoPointCalculator::weightInput
(
    const point p,
    const point& centre,
    const scalar R
) const
{
    point projected = centre + R*(p - centre) / ((p & centre) + EPSILON);
    return mag(p - projected);
}

}} // End namespace Foam::FrontTracking 

// ************************************************************************* //
