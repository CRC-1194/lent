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
    A class for for the RBF interpolation over a mesh.  

SourceFiles
    rbfCellsInterpolationEigen.C

\*---------------------------------------------------------------------------*/

#ifndef rbfCellsInterpolationEigen_H
#define rbfCellsInterpolationEigen_H

#include "rbfConfig.H"
#include "rbfInterpolationEigen.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBF
{

// FIXME: Move to a separate file. TM.
// FIXME: Optimize for construct/copy/move/resize operations.
struct cellRbfStencil
{
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

    label cellLabel_;
    std::vector<label>  pointLabels_; 
    std::vector<label>  cellNeiLabels_;
    int size_; 
};

/*---------------------------------------------------------------------------*\
                         Class rbfCellsInterpolationEigen Declaration
\*---------------------------------------------------------------------------*/

// FIXME: 
// - Check if the mesh is equidistant and in this case, provide only a single 
//   factorized matrix. 
// - Read the stencilType from a config dictionary.  
class rbfCellsInterpolationEigen final
{

    constexpr static int defaultStencilSize_ = 15; // BCC_FVM stencil size for a hex cell. 
public:

    // Constructors
    rbfCellsInterpolationEigen() = default;

    // FIXME: Dict constructor that reads stencilType and support type.
    template<typename Kernel>
    rbfCellsInterpolationEigen
    (
        const fvMesh& mesh, 
        const Kernel& kernel, 
        const stencilType stencil = stencilType::BCC,  
        const supportType support = supportType::GLOBAL_SUPPORT
    ); 

    // Member Functions
    void calcStencils 
    (
        const fvMesh& mesh,
        const stencilType stencil = stencilType::BCC
    );
    
    void factorize();

    void interpolate
    (
        const volScalarField& vf, 
        const pointScalarField& pf
    );
    

    // TODO: Inline and measure effect.
    template<typename Point>
    double evaluate
    (
        Point const& evalPoint, 
        const label cellI
    );

    // TODO: Inline and measure effect.
    void setCellRbfValues
    (
        const volScalarField& vf, 
        const pointScalarField& pf, 
        const label cellI 
    );

    
private:

    // Cell RBF stencils.
    std::vector<cellRbfStencil> cellRbfStencils_; 
    
    // Cell-centered RBF point lists.
    std::vector<std::vector<Foam::point>> cellRbfPoints_; 

    // Cell-centered RBF value lists.
    std::vector<Eigen::VectorXd> cellRbfValues_; 
    
    // RBF interpolants.
    std::vector<rbfInterpolationEigen> cellRbfs_; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBF 
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rbfCellsInterpolationEigenTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
