/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
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

#include "rbfCellsInterpolationEigen.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBF
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void rbfCellsInterpolationEigen::calcStencils
(
    const fvMesh& mesh, 
    const stencilType stencil 
)
{
    // FIXME: emplace_back will not work with re-calculation, fix this. 
    // Resizing to zero, while maintaining same capacity as before when 
    // re-executing. TM
    
    // Append cell-centered data for both BCC and BCC-FVM stencil. 
    const auto& cellCenters = mesh.C(); 
    forAll (cellCenters, cellI)
    {
        cellRbfStencils_[cellI].setCellLabel(cellI);
        cellRbfPoints_[cellI].emplace_back(cellCenters[cellI]);
    }

    // Add point-data to for both the BCC and BCC-FVM stencil.
    const auto& pointCells = mesh.pointCells(); 
    const auto& meshPoints = mesh.points(); 
    forAll (pointCells, pointI)
    {
        const auto& pointCellLabels = pointCells[pointI]; 
        forAll(pointCellLabels, cellI)
        {
            const auto cellLabel = pointCellLabels[cellI]; 
            cellRbfPoints_[cellLabel].emplace_back(meshPoints[pointI]);
            cellRbfStencils_[cellLabel].addPointLabel(pointI);
        }
    }
    // TODO: Parallel implementation for point data. TM.

    if (stencil == stencilType::BCC_FVM)
    {
        // Add face-neighbor data to the BCC-FVM cell stencil.
        const auto& own = mesh.owner(); 
        const auto& nei = mesh.neighbour();
        forAll(own, faceI)
        {
            const auto ownCell = own[faceI]; 
            const auto neiCell = nei[faceI];  
            cellRbfPoints_[ownCell].emplace_back(cellCenters[neiCell]);
            cellRbfPoints_[neiCell].emplace_back(cellCenters[ownCell]);

            cellRbfStencils_[ownCell].addNeiLabel(neiCell); 
            cellRbfStencils_[neiCell].addNeiLabel(ownCell); 
        }
        // TODO: Parallel implementation for processor fv patches. TM.
    }
}

void rbfCellsInterpolationEigen::factorize()
{
    for (decltype(cellRbfPoints_.size()) cellI = 0; cellI < cellRbfPoints_.size(); ++cellI)
        cellRbfs_[cellI].factorize(cellRbfPoints_[cellI]);
}

void rbfCellsInterpolationEigen::setCellRbfValues
(
    const volScalarField& vf, 
    const pointScalarField& pf, 
    const label cellI 
)
{
    const cellRbfStencil& rbfStencil = cellRbfStencils_[cellI]; 
    Eigen::VectorXd& rbfValues_ = cellRbfValues_[cellI];  
    rbfValues_.resize(rbfStencil.size()); 

    // 0 place is used for the cell value given by the RBF cell label.
    rbfValues_[0] = vf[rbfStencil.cellLabel()];  

    // [1, NcellPoints] are the RBF cell corner point labels. 
    const auto& rbfPointLabels = rbfStencil.pointLabels();  
    using ptSizeType = decltype(rbfPointLabels.size()); 
    for (ptSizeType pI = 0; pI < rbfPointLabels.size(); ++pI)
        rbfValues_[pI + 1] = pf[rbfPointLabels[pI]]; 

    //const auto& rbfCellNeiLabels = rbfStencil.cellNeiLabels(); 
    const auto rbfCellNeiLabels = rbfStencil.cellNeiLabels(); 
    using clSizeType = decltype(rbfCellNeiLabels.size()); 
    const auto shift = 1 + rbfPointLabels.size(); 
    for (clSizeType pI = 0; pI < rbfCellNeiLabels.size(); ++pI)
        rbfValues_[pI + shift] = pf[rbfPointLabels[pI]]; 
}

void rbfCellsInterpolationEigen::interpolate
(
    const volScalarField& vField,
    const pointScalarField& pField
)
{
    
    // For each cell 
    forAll(vField, cellI)
    {
        setCellRbfValues(vField, pField, cellI);
        cellRbfs_[cellI].interpolate(cellRbfPoints_[cellI], cellRbfValues_[cellI]); 
    }
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBF 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
