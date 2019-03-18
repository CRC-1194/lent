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

Author: 
    Tomislav Maric, maric@mma.tu-darmstadt.de
    TU Darmstadt, Mathematical modeling and analysis group 

\*---------------------------------------------------------------------------*/

#include "rbfInterpolationEigen.H"
#include "rbfGeometry.H"
#include "vector.H"

namespace Foam { namespace RBF {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Kernel>
rbfInterpolationEigen<Kernel>::rbfInterpolationEigen()
:
    kernel_(), 
    support_(supportType::GLOBAL_SUPPORT), 
    supportRadius_(1.0),
    matrix_(), 
    matrixFactorization_(), 
    source_(), 
    rbfCoeffs_()
{}

template<typename Kernel> 
rbfInterpolationEigen<Kernel>::rbfInterpolationEigen(const supportType support)
:
    kernel_(), 
    support_(support), 
    supportRadius_(1.0),
    matrix_(), 
    matrixFactorization_(), 
    source_(), 
    rbfCoeffs_()
{}

template<typename Kernel> 
template<typename Points> 
rbfInterpolationEigen<Kernel>::rbfInterpolationEigen
(
    Points const& points, 
    const supportType support 
)
:
    kernel_(), 
    support_(support), 
    supportRadius_(1.0),
    matrix_(points.size() + Npoly_, points.size() + Npoly_), 
    matrixFactorization_(), 
    source_(Npoly_ + points.size()),
    rbfCoeffs_(Npoly_ + points.size())
{
    factorize(points); 
}

template<typename Kernel>
template<typename Points, typename Values> 
rbfInterpolationEigen<Kernel>::rbfInterpolationEigen
(
    Points const& points, 
    Values const& values, 
    const supportType support 
)
:
    kernel_(), 
    support_(support), 
    supportRadius_(1.0),
    matrix_(points.size() + Npoly_, points.size() + Npoly_), 
    matrixFactorization_(), 
    source_(Npoly_ + points.size()),
    rbfCoeffs_(Npoly_ + points.size())
{
    factorize(points); 
    solve(points, values);
}

// * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * * //

template<typename Kernel>
template<typename Points> 
void rbfInterpolationEigen<Kernel>::factorize(Points const& points)
{
    matrix_.resize(points.size() + Npoly_, points.size() + Npoly_); 

    supportRadius_ = 1.0;

    if (support_ == supportType::MAX_CENTROID_RADIUS)
        supportRadius_ = Geometry::maxCentroidRadius(points);  

    if (support_ == supportType::MEAN_CENTROID_RADIUS)
        supportRadius_ = Geometry::meanCentroidRadius(points);  

    if (support_ == supportType::RMS_RADIUS)
        supportRadius_ = Geometry::rmsRadius(points);  

    const auto Npoints = points.size(); 

    using sizeType = decltype(points.size()); 

    // Fill the RBF kernel coefficients block.
    for (sizeType j = 0; j < Npoints; ++j)
        for (sizeType i = 0; i < Npoints ; ++i)
        {
            const double dist = Geometry::distance(points[i], points[j]); 
            matrix_(i,j) = kernel_(dist, supportRadius_);
        }

    // Fill the linear polynomial coefficients block.
    auto P = matrix_.topRightCorner(Npoints, Npoly_);  
    auto Pt = matrix_.bottomLeftCorner(Npoly_, Npoints);
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
    matrix_.bottomRightCorner(Npoly_, Npoly_) = MatrixType::Zero(Npoly_, Npoly_); 

    // Factorize the matrix.
    matrixFactorization_.compute(matrix_); 
}

template<typename Kernel>
template<typename Points, typename Values>
void rbfInterpolationEigen<Kernel>::solve
(
    Points const& points, 
    Values const& values
) 
{
    // Extend the "values" source term with 0s of the linear polynomial. 
    source_.resize(Npoly_ + points.size()); 
    source_.head(points.size()) = values; 
    source_.tail(Npoly_) = Eigen::VectorXd::Zero(Npoly_); 

    // Solve for the RBF coefficients.
    rbfCoeffs_ = matrixFactorization_.solve(source_);
}

template<typename Kernel> 
template<typename Point, typename Points> 
double rbfInterpolationEigen<Kernel>::value(
    Point const& point, 
    Points const& points
)  const
{
    double result = 0.; 

    const auto Npoints = points.size(); 

    using sizeType = decltype(points.size()); 
    for (sizeType i = 0; i < Npoints; ++i)
    {
        const double dist = Geometry::distance(point, points[i]); 
        result += rbfCoeffs_[i] * 
                  kernel_(dist, supportRadius_); 
    }

    result += rbfCoeffs_[Npoints]; 

    for (char i = 1; i < Npoly_; ++i)
        result += rbfCoeffs_[i + Npoints] * point[i - 1]; 

    return result;
}

template<typename Kernel> 
template<typename Vector, typename Points> 
Vector rbfInterpolationEigen<Kernel>::grad(
    Vector const& point, 
    Points const& points
)  const
{
    Vector result (0,0,0);  

    const auto Npoints = points.size(); 

    using sizeType = decltype(points.size()); 
    for (sizeType i = 0; i < Npoints; ++i)
        result += rbfCoeffs_[i] * kernel_.template grad<Vector>(point, points[i]); 

    result += Vector(
        rbfCoeffs_[Npoints + 1], 
        rbfCoeffs_[Npoints + 2],  
        rbfCoeffs_[Npoints + 3]  
    );

    return result;
}

}} // End namespace Foam::RBF 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

