/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Tomislav Maric, TU Darmstadt
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

#ifndef rbfIsoInterpolationTemplates_H
#define rbfIsoInterpolationTemplates_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace RBF {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Kernel>
rbfIsoInterpolation::rbfIsoInterpolation
(
    const fvMesh& mesh, 
    const Kernel& kernel,
    const stencilType stencil,
    const supportType support
)
    :
        cellRbfStencils_(mesh.nCells()),
        cellRbfPoints_(mesh.nCells()),
        cellRbfValues_(mesh.nCells()),
        cellRbfs_(mesh.nCells(), 
                  rbfInterpolation<Kernel>())

{
    calcStencils(mesh, stencil); 
    factorize(kernel, support); 
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<typename Kernel>
void rbfIsoInterpolation::factorize
(
    Kernel const& kernel, 
    const supportType support
)
{
    for (decltype(cellRbfPoints_.size()) cellI = 0; cellI < cellRbfPoints_.size(); ++cellI)
        cellRbfs_[cellI].factorize(cellRbfPoints_[cellI], kernel, support);
}

template<typename Point> 
double rbfIsoInterpolation::evaluate
(   
    Point const& evalPoint, 
    const label cellI
)
{
    return cellRbfs_[cellI].evaluate
    (
        evalPoint,
        cellRbfPoints_[cellI], 
        cellRbfValues_[cellI]
    ); 
}

}} // End namespace Foam::RBF 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
