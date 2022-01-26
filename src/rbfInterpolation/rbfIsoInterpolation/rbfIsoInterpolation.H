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

#ifndef rbfIsoInterpolation_H
#define rbfIsoInterpolation_H

#include "rbfConfig.H"
#include "rbfInterpolationEigen.H"
#include "rbfCellsInterpolationEigen.H"
#include "fvMesh.H"
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace RBF {

/*---------------------------------------------------------------------------*\
                         Class rbfIsoInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<typename Kernel>
class rbfIsoInterpolation final
{

private:

    // Cell RBF stencils.
    std::vector<cellRbfStencil> cellRbfStencils_; 
    
    // Cell-centered RBF points.
    std::vector<std::vector<Foam::point>> cellRbfPoints_; 

    // Cell-centered RBF values. 
    std::vector<::VectorXd> cellRbfValues_; 
    
    // Cell-centered RBF interpolants.
    std::vector<rbfInterpolation<Kernel>> cellRbfs_; 

public:

    // Constructors
    rbfIsoInterpolation() = default;

    rbfIsoInterpolation
    (
        const volScalarField& vField, 
        const pointScalarField& pField, 
        const stencilType stencil = stencilType::BCCC,  
        const supportType support = supportType::GLOBAL_SUPPORT
    ); 

    // Member Functions
    template<typename CellLabels>
    void identifyIsoCells(
        CellLabels& cellLabels,  
        const volScalarField& vField, 
        const pointScalarField& pField, 
        const scalar isoValue
    );

    void calcStencils 
    (
        const volScalarField& vField, 
        const pointScalarField& pField,
        const stencilType stencil
    );
    
    /*
    void factorize(); 

    void solve 
    (
        const volScalarField& vf, 
        const pointScalarField& pf
    );
    
    template<typename Point>
    double value (Point const& evalPoint, const label cellI) const;

    template<typename Point>
    double value (const label cellI, Point const& evalPoint) const;

    template<typename Vector, typename Point> 
    Vector grad(const label cellI, Point const& point) const; 

    template<typename Vector, typename Point> 
    Vector grad(Point const& point, const label cellI) const; 

    void setCellRbfValues
    (
        const volScalarField& vf, 
        const pointScalarField& pf, 
        const label cellI 
    );

    const rbfInterpolation<Kernel>& rbf(std::size_t cellI) const
    {
        return cellRbfs_[cellI]; 
    }

    const std::vector<Foam::point>& rbfPoints(std::size_t cellI) const
    {
        return cellRbfPoints_[cellI]; 
    }

    std::string rbfName() const
    {
        return cellRbfs_.front().name(); 
    }
    */

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::RBF

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rbfIsoInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
