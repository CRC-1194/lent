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
    Foam::isoPointCalculator

Description
    Calculates points on an iso-surface \Gamma := {x, f(x) = s}, where 
    s is the isovalue using an RBF interpolation. 

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    isoPointCalculatorI.H
    isoPointCalculator.C

\*--------------------------------------------------------------------------*/

#include "isoPointCalculator.H"
#include "volMesh.H"
#include <cassert> 

namespace Foam { namespace FrontTracking {

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool isoPointCalculator::isoValueInInterval
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

bool isoPointCalculator::equalUnderTolerance(
    const scalar s1, 
    const scalar s2, 
    const scalar tol
) const
{
    return (mag(s1 - s2) < tol * max(max(mag(s1), mag(s2)), 1.0));
}

void isoPointCalculator::calcEdgePoints
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
    edgeLabels_.resize(mesh.nEdges());
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
            edgeLabels_[edgeI] = edgePoints_.size() - 1; 
        }
    }
}

void isoPointCalculator::calcCellPoints
(
    const volScalarField& cellPhi, 
    const pointScalarField& pointPhi, 
    const scalar isoValue 
)
{
    calcEdgePoints(pointPhi, isoValue); 
    
    const auto& mesh = cellPhi.mesh();
    const auto& meshCellEdgesLists = mesh.cellEdges(); 

    // FIXME: Use reserve in the constructor. TM.
    cellPoints_.resize(0);
    cellLabels_.resize(mesh.nCells());
    // FIXME: Can also be a list of booleans. TM.
    std::fill(cellLabels_.begin(), cellLabels_.end(), -1);

    forAll(meshCellEdgesLists, cellI)
    {
        const auto& cellEdgeLabels = meshCellEdgesLists[cellI]; 
        label nEdgePoints = 0;
        auto cellPoint = vector(0,0,0);
        forAll(cellEdgeLabels, edgeI)  
        {
            const label edgeGl = cellEdgeLabels[edgeI]; 
            if (edgeLabels_[edgeGl] > 0)
            {
                cellPoint += edgePoints_[edgeLabels_[edgeGl]];  
                ++nEdgePoints; 
            }
        }
        if (nEdgePoints > 0)
        {
            cellPoints_.push_back(cellPoint / nEdgePoints); 
            cellLabels_[cellI] = cellI;  
        }
    }
}

}} // End namespace Foam::FrontTracking 

// ************************************************************************* //
