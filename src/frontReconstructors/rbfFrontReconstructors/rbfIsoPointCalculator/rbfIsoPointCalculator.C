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
    Foam::rbfIsoPointCalculator

Description
    Calculates points on an iso-surface \Gamma := {x, f(x) = s}, where 
    s is the isovalue using an RBF interpolation. 

Authors
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis Group, TU Darmstadt

SourceFiles
    rbfIsoPointCalculatorI.H
    rbfIsoPointCalculator.C

\*--------------------------------------------------------------------------*/

#include "rbfIsoPointCalculator.H"
#include "volMesh.H"
#include <cassert> 

namespace Foam { namespace FrontTracking {

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<typename MeshRbfs>
void rbfIsoPointCalculator::rbfCorrectContourPoint(
    point& contourPoint, 
    label pointI,
    MeshRbfs const& meshRbfs 
)
{
    // Compute the value of the signed distance RBF at point pointI.
    double rbfVal = meshRbfs.value(pointCellLabels_[pointI], contourPoint);
        vector rbfGrad = meshRbfs.template grad<vector>(
            pointCellLabels_[pointI], 
            contourPoint
        );

    //contourPointRbfVals_.push_back(rbfVal); 
    //contourPointRbfGrads_.push_back(rbfGrad); 

    //if (rbfVal > 1e-3) // Propagate isoValue up to here? TM.
    //{
    contourPoint += -(rbfVal * rbfGrad) / (rbfGrad & rbfGrad);
    //}

    // BEGIN Debugging, Sphere test center=(0.,0.,0.), radius = 0.5. 
    // Test application: lentTestRbfReconstruction.
    //vector exactGrad(
        //2 * (contourPoint[0]), 
        //2 * (contourPoint[1]),
        //2 * (contourPoint[2])
    //);
    //scalar exactValue(Foam::sqr(contourPoint[0]) + Foam::sqr(contourPoint[1]) + Foam::sqr(contourPoint[2]) -  0.25);
    //const scalar exactGradMagSqr = exactGrad & exactGrad; 
    //contourPoint += -(exactValue * exactGrad) / exactGradMagSqr;
    // END Debugging
}

template<typename MeshRbfs>
void rbfIsoPointCalculator::rbfCorrectContourPoints(MeshRbfs const& meshRbfs)
{
    for(decltype(contourPoints_.size()) pointI = 0; pointI < contourPoints_.size(); ++pointI)
        rbfCorrectContourPoint(contourPoints_[pointI], static_cast<Foam::label>(pointI), meshRbfs);
}

template<typename MeshRbfs>
void rbfIsoPointCalculator::calcContourPoints
(
    const volScalarField& cellPhi, 
    const pointScalarField& pointPhi,
    MeshRbfs const&, 
    const scalar isoValue 
)
{
    calcEdgePoints(pointPhi, isoValue);
    
    const auto& mesh = cellPhi.mesh();
    const auto& meshCellEdgesLists = mesh.cellEdges(); 

    contourPoints_.resize(0);
    pointCellLabels_.resize(0);
    //contourPointRbfVals_.resize(0); 
    //contourPointRbfGrads_.resize(0);
    //cellLabels_.resize(mesh.nCells());
    //std::fill(cellLabels_.begin(), cellLabels_.end(), -1);

    // Compute dual-contouring points using averaging and RBF-based projection
    // onto the RBF implicit surface.
    forAll(meshCellEdgesLists, cellI)
    {
        const auto& cellEdgeLabels = meshCellEdgesLists[cellI]; 
        label nEdgePoints = 0;
        auto contourPoint = vector(0,0,0);
        forAll(cellEdgeLabels, edgeI)  
        {
            const auto edgeGl = static_cast<unsigned long>(cellEdgeLabels[edgeI]); 
            if (edgeLabels_[edgeGl] > 0)
            {
                contourPoint += edgePoints_[static_cast<unsigned long>(edgeLabels_[edgeGl])];  
                ++nEdgePoints; 
            }
        }
        if (nEdgePoints > 0)
        {
            contourPoint /= nEdgePoints;  
            // This is the initial contourPoint. It needs to be corrected
            // using RBF interpolation. 
            // FIXME: Enable this.
            //rbfCorrectContourPoint(contourPoint, cellI, meshRbfs); 
           
            contourPoints_.push_back(contourPoint); 
            pointCellLabels_.push_back(cellI);
            //cellLabels_[cellI] = cellI;  
        }
    }
}


template<typename MeshRbfs>
std::vector<double> rbfIsoPointCalculator::contourPointValues(MeshRbfs const& meshRbfs) const
{
    std::vector<double> values(contourPoints_.size());
    for(decltype(contourPoints_.size()) pointI = 0; pointI < contourPoints_.size(); ++pointI)
        values[pointI] = meshRbfs.value(pointCellLabels_[pointI], contourPoints_[pointI]);

    return values;
}

}} // End namespace Foam::FrontTracking 

// ************************************************************************* //
