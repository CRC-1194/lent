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
    Foam::rbfCellsInterpolationEigen

Description
    A class for for the RBF interpolation that is cell-centered. 
    It relies on stencils centered in the cell: 

    * BCC - body centered cubic 
    * BCCC - body centered cubic + first face neighbors 

SourceFiles
    rbfCellsInterpolationEigen.C

\*---------------------------------------------------------------------------*/

#ifndef rbfCellsInterpolationEigen_H
#define rbfCellsInterpolationEigen_H

#include "rbfConfig.H"
#include "rbfInterpolationEigen.H"
#include "fvMesh.H"
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace RBF {

/*---------------------------------------------------------------------------*\
                         Class cellRbfStencil Definition 
\*---------------------------------------------------------------------------*/

class cellRbfStencil
{
    public: 

        cellRbfStencil()
            : 
                cellLabel_(-1), 
                pointLabels_(), 
                cellNeiLabels_(),
                size_(0)
        {}

        void setCellLabel(const label cellI)
        {
            cellLabel_ = cellI; 
            ++size_;
        }

        void addPointLabel(const label pointI) 
        {
            pointLabels_.push_back(pointI); 
            ++size_; 
        }

        void addNeiLabel(const label cellI) 
        {
            cellNeiLabels_.push_back(cellI); 
            ++size_; 
        }

        int size() const
        {
            return size_;  
        }

        label cellLabel() const
        {
            return cellLabel_; 
        }

        const std::vector<label>& pointLabels() const
        {
            return pointLabels_; 
        }

        const std::vector<label>& cellNeiLabels() const
        {
            return cellNeiLabels_;
        }

    private:

        label cellLabel_;
        std::vector<label>  pointLabels_; 
        std::vector<label>  cellNeiLabels_;
        int size_; 
};

/*---------------------------------------------------------------------------*\
                         Class rbfCellsInterpolationEigen Declaration
\*---------------------------------------------------------------------------*/

template<typename Kernel>
class rbfCellsInterpolationEigen final
{

private:

    // Cell RBF stencils.
    std::vector<cellRbfStencil> cellRbfStencils_; 
    
    // Cell-centered RBF points.
    std::vector<std::vector<Foam::point>> cellRbfPoints_; 

    // Cell-centered RBF values. 
    std::vector<Eigen::VectorXd> cellRbfValues_; 
    
    // Cell-centered RBF interpolants.
    std::vector<rbfInterpolationEigen<Kernel>> cellRbfs_; 

public:

    // Constructors
    rbfCellsInterpolationEigen() = default;

    rbfCellsInterpolationEigen
    (
        const fvMesh& mesh, 
        const stencilType stencil = stencilType::BCCC,  
        const supportType support = supportType::GLOBAL_SUPPORT
    ); 

    // Member Functions
    void calcStencils 
    (
        const fvMesh& mesh, 
        const stencilType stencil 
    );
    
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

    const rbfInterpolationEigen<Kernel>& rbf(std::size_t cellI) const
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

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::RBF

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rbfCellsInterpolationEigen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
