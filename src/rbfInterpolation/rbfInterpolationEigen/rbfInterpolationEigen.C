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
    Foam::rbfInterpolationEigen

Description
    RBF interpolation extended with a linar polynomial. 

    Linear system solution implemented using the Eigen library.  

    The interpolant is constructed using a colllection of points and a
    collection of values. 

    The interpolation uses the Eigen library is to find the solution to the 
    linear system. LDLT decomposition is used for the system solution. 

    It is not applicable for a large number of points/values. 

    NO LOGICAL CHECKS ARE PERFORMED FOR EFFICIENCY REASONS

        - To evaluate the interpolant at a point, the RBF system needs to
          be first factorized, then interpolated (solved).

        - This allows the user to first factorize the system, and then 
          re-use the factorization with different source terms, and then
          evaluate the interpolant at a point in the stencil.

        1. It is assumed that the factorization, interpolation and evaluation
           is performed with the same set of nodal points. 
        2. It is assumed that the nodal points and nodal values have the same
           length. 
        3. It is assumed that the system has bee factorized and solved for the
           evaluation execution. 
        4. It is assumed that the evaluation point is "within" the stencil: 
           evaluation far away from the nodal points results in a higher 
           interpolation error.

SourceFiles
    rbfInterpolationEigen.H
    rbfInterpolationEigen.C

\*---------------------------------------------------------------------------*/

#include "rbfInterpolationEigen.H"
#include "rbfGeometry.H"

namespace Foam
{
namespace RBF 
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rbfInterpolationEigen::rbfInterpolationEigen()
:
    A_(), 
    matrixFactorization_(), 
    rbfCoeffs_(), 
    kernel_(),
    support_(supportType::GLOBAL_SUPPORT), 
    supportRadius_(1.0)
{}

template<typename Kernel> 
rbfInterpolationEigen::rbfInterpolationEigen
(
    const Kernel& kernel, 
    supportType support
)
:
    A_(), 
    matrixFactorization_(), 
    rbfCoeffs_(), 
    kernel_(kernel), 
    support_(support),
    supportRadius_(1.0)
{}

template<typename Points, typename Kernel> 
rbfInterpolationEigen::rbfInterpolationEigen
(
    Points const& points, 
    const Kernel& kernel, 
    const supportType support 
)
:
    A_(points.size() + Npoly_, points.size() + Npoly_), 
    matrixFactorization_(), 
    rbfCoeffs_(), 
    kernel_(kernel), 
    support_(support),
    supportRadius_(1.0)
{
    factorize(points); 
}

template<typename Points, typename Values, typename Kernel> 
rbfInterpolationEigen::rbfInterpolationEigen
(
    Points const& points, 
    Values const& values, 
    const Kernel& kernel, 
    const supportType support 
)
:
    A_(points.size() + Npoly_, points.size() + Npoly_), 
    matrixFactorization_(), 
    rbfCoeffs_(), 
    kernel_(kernel), 
    support_(support),
    supportRadius_(1.0)
{
    factorize(points); 
    interpolate(points, values);
}

// * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

template<typename Points> 
void rbfInterpolationEigen::factorize(Points const& points)
{
    A_.resize(points.size() + Npoly_, points.size() + Npoly_); 

    if (support_ == supportType::MAX_CENTROID_RADIUS)
        supportRadius_ = Geometry::maxCentroidRadius(points);  

    if (support_ == supportType::MEAN_CENTROID_RADIUS)
        supportRadius_ = Geometry::meanCentroidRadius(points);  

    if (support_ == supportType::RMS_RADIUS)
        supportRadius_ = Geometry::rmsRadius(points);  

    const auto Npoints = points.size(); 

    // Fill the RBF kernel coefficients block.
    for (sizeType j = 0; j < Npoints; ++j)
        for (sizeType i = 0; i < Npoints ; ++i)
            A_(i,j) = kernel_(Geometry::distance(points[i], points[j]), supportRadius_);

    // Fill the linear polynomial coefficients block.
    auto P = A_.topRightCorner(Npoints, Npoly_);  
    auto Pt = A_.bottomLeftCorner(Npoly_, Npoints);
    for(sizeType i = 0; i < Npoints; ++i)
    {
        P(i, 0) = 1.; 
        Pt(0, i) = 1.;
        for(char k = 1; k < Npoly_; ++k)
        {
            P(i,k) = points[i][k-1]; 
            Pt(k,i) = points[i][k-1]; 
        }
    }
    // Fill the lower right block with 0s.
    A_.bottomRightCorner(Npoly_, Npoly_) = MatrixType::Zero(Npoly_, Npoly_); 

    // Assemble the matrix and factorize it in place. 
    //matrixFactorization_ = A_.llt(); 
    //matrixFactorization_ = A_.ldlt(); 
    //matrixFactorization_ = A_.householderQr(); 
    matrixFactorization_ = A_.partialPivLu(); 
}

template<typename Points, typename Kernel> 
void rbfInterpolationEigen::factorize(Points const& points, const Kernel& kernel)
{
    kernel_ = kernel;   
    factorize(points); 
}

template<typename Points, typename Values>
void rbfInterpolationEigen::interpolate
(
    Points const& points, 
    Values const& values
) 
{
    // Extend the "values" source term with 0s of the linear polynomial. 
    Eigen::VectorXd source(Npoly_ + points.size()); 
    source.head(points.size()) = values; 
    source.tail(Npoly_) = Eigen::VectorXd::Zero(Npoly_); 

    // Solve for the RBF coefficients.
    rbfCoeffs_ = matrixFactorization_.solve(source);
}

template<typename Points, typename Values, typename Kernel>
void rbfInterpolationEigen::interpolate
(
    Points const& points, 
    Values const& values, 
    const Kernel& kernel
)
{
    kernel_ = kernel; 
    interpolate(points, values);
}

template<typename Point, typename Points, typename Values> 
double rbfInterpolationEigen::evaluate(
    Point const& point, 
    Points const& points, 
    Values const& values
) 
{
    double result = 0.; 

    const auto Npoints = points.size(); 

    for (sizeType i = 0; i < Npoints; ++i)
        result += rbfCoeffs_[i] * 
                  kernel_(Geometry::distance(point, points[i]), supportRadius_); 

    result += rbfCoeffs_[Npoints]; 

    for (char i = 1; i < Npoly_; ++i)
        result += rbfCoeffs_[i + Npoints] * point[i - 1]; 

    return result;
}


} // End namespace RBF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

