/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, TU Darmstadt
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

Class
    Foam::rbfIsoInterpolation

Description
    A class for for the RBF interpolation that is cell-centered. 
    It relies on stencils centered in the cell: 

    * BCC - body centered cubic 
    * BCCC - body centered cubic + first face neighbors 

    The rbfIsoInterpolation constructs the RBF interpolation only in those cells
    that are intersected by an iso-contour. Use rbfCellsInterpolationEigen when 
    an RBF interpolation is required in each cell.

SourceFiles
    rbfIsoInterpolation.H
    rbfIsoInterpolation.C
    rbfIsoInterpolationTemplates.H

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace RBF {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Kernel>
rbfIsoInterpolation<Kernel>::rbfIsoInterpolation
(
    const volScalarField& vfield, 
    const pointScalarField& pfield, 
    const stencilType stencil,
    const supportType support
)
    :
        isoValue_(isoValue),
        cellRbfStencils_(),
        cellRbfPoints_(),
        cellRbfValues_(),
        cellRbfs_() 
{
    calcStencils(vField, pField, isoValue, stencil); 
    factorize(); 
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<typename Kernel>
void rbfIsoInterpolation<Kernel>::calcStencils
(
    const volScalarField& vField, 
    const pointScalarField& pField,
    const scalar isoValue, 
    const stencilType stencil
)
{
    /*
    // FIXME: emplace_back will not work with re-calculation! TM. 
    
    // Append cell-centered data for both BCC and BCCC stencil. 
    const auto& cellCenters = mesh.C(); 
    forAll (cellCenters, cellI)
    {
        cellRbfStencils_[cellI].setCellLabel(cellI);
        cellRbfPoints_[cellI].emplace_back(cellCenters[cellI]);
    }

    // Add point-data to for both the BCC and BCCC stencil.
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

    if (stencil == stencilType::BCCC)
    {
        // Add face-neighbor data to the BCCC cell stencil.
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
    */
}

template<typename Kernel>
void rbfIsoInterpolation<Kernel>::setCellRbfValues
(
    const volScalarField& vf, 
    const pointScalarField& pf, 
    const label cellI 
)
{
    const cellRbfStencil& rbfStencil = cellRbfStencils_[cellI]; 
    ::VectorXd& rbfValues_ = cellRbfValues_[cellI];  
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
        rbfValues_[pI + shift] = vf[rbfCellNeiLabels[pI]]; 
}

template<typename Kernel>
void rbfIsoInterpolation<Kernel>::factorize()
{
    using sizeType = decltype(cellRbfPoints_.size());
    for (sizeType cellI = 0; cellI < cellRbfPoints_.size(); ++cellI)
        cellRbfs_[cellI].factorize(cellRbfPoints_[cellI]);
}

template<typename Kernel>
template<typename Point> 
double rbfIsoInterpolation<Kernel>::value
(   
    Point const& evalPoint, 
    const label cellI
) const
{
    return cellRbfs_[cellI].value(evalPoint, cellRbfPoints_[cellI]); 
}

template<typename Kernel>
template<typename Point> 
double rbfIsoInterpolation<Kernel>::value
(   
    const label cellI,
    Point const& evalPoint 
) const
{
    return this->value(evalPoint, cellI);
}

template<typename Kernel>
template<typename Vector, typename Point> 
Vector rbfIsoInterpolation<Kernel>::grad
(   
    const label cellI,
    Point const& evalPoint
) const
{
    return cellRbfs_[cellI].template grad<Vector>(
        evalPoint, 
        cellRbfPoints_[cellI]
    ); 
}

template<typename Kernel>
template<typename Vector, typename Point> 
Vector rbfIsoInterpolation<Kernel>::grad
(   
    Point const& evalPoint,
    const label cellI
) const
{
    return this->template grad<Vector>(cellI, evalPoint);
}

template<typename Kernel>
void rbfIsoInterpolation<Kernel>::solve
(
    const volScalarField& vField,
    const pointScalarField& pField
)
{
    forAll(vField, cellI)
    {
        setCellRbfValues(vField, pField, cellI);
        cellRbfs_[cellI].solve(cellRbfPoints_[cellI], cellRbfValues_[cellI]); 
    }
}

}} // End namespace Foam::RBF

// ************************************************************************* //
