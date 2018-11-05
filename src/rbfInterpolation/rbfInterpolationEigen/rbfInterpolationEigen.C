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

SourceFiles
    rbfInterpolationEigen.H
    rbfInterpolationEigen.C

\*---------------------------------------------------------------------------*/

#include "rbfInterpolationEigen.H"
#include "rbfGeometry.H"

#include <limits>

namespace RBF 
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rbfInterpolationEigen::rbfInterpolationEigen()
:
    matrixFactorization_(), 
    rbfCoeffs_(), 
    kernel_(),
    supportRadius_(1.0),
    factorized_(false), 
    interpolated_(false)
{}

template<typename Points, typename Kernel> 
rbfInterpolationEigen::rbfInterpolationEigen
(
    Points const& points, 
    Kernel kernel, 
    supportType support 
)
:
    matrixFactorization_(), 
    rbfCoeffs_(), 
    kernel_(kernel), 
    supportRadius_(1.0), 
    factorized_(false),
    interpolated_(false)
{
    factorize(points, support); 
}

template<typename Points, typename Values, typename Kernel> 
rbfInterpolationEigen::rbfInterpolationEigen
(
    Points const& points, 
    Values const& values, 
    Kernel kernel, 
    supportType support 
)
:
    matrixFactorization_(), 
    rbfCoeffs_(), 
    kernel_(kernel), 
    supportRadius_(1.0),
    factorized_(false),
    interpolated_(false)
{
    interpolate(points, values, support);
}

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //
//
template<typename Points> 
Eigen::MatrixXd rbfInterpolationEigen::matrix
(
    Points const& points, 
    supportType support 
)
{
    // NRVO initialization.
    Eigen::MatrixXd A (points.size() + Npoly_, points.size() + Npoly_);

    if (support == supportType::MAX_CENTROID_RADIUS)
        supportRadius_ = Numeric::maxCentroidRadius(points);  

    if (support == supportType::MEAN_CENTROID_RADIUS)
        supportRadius_ = Numeric::meanCentroidRadius(points);  

    if (support == supportType::RMS_RADIUS)
        supportRadius_ = Numeric::rmsRadius(points);  

    const auto Npoints = points.size(); 

    // Fill the RBF kernel coefficients block.
    for (sizeType j = 0; j < Npoints; ++j)
        for (sizeType i = 0; i < Npoints ; ++i)
            A(i,j) = kernel_(Numeric::distance(points[i], points[j]), supportRadius_);

    // Fill the linear polynomial coefficients block.
    auto P = A.topRightCorner(Npoints, Npoly_);  
    auto Pt = A.bottomLeftCorner(Npoly_, Npoints);
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
    A.bottomRightCorner(Npoly_, Npoly_) = Eigen::MatrixXd::Zero(Npoly_, Npoly_); 
    
    return A; 
}

template<typename Points> 
void rbfInterpolationEigen::factorize
(
    Points const& points, 
    supportType support 
)
{
    Eigen::MatrixXd A = matrix(points, support);  

    matrixFactorization_ = A.colPivHouseholderQr(); 

    factorized_ = true;
}

template<typename Points, typename Values>
void rbfInterpolationEigen::interpolate
(
    Points const& points, 
    Values const& values, 
    supportType support 
) 
{
    if (! factorized_)
        factorize(points, support);

    // Extend the "values" source term with 0s of the linear polynomial. 
  Eigen::VectorXd source(Npoly_ + points.size()); 
    source.head(points.size()) = values; 
    source.tail(Npoly_) = Eigen::VectorXd::Zero(Npoly_); 

    // Solve for the RBF coefficients.
    rbfCoeffs_ = matrixFactorization_.solve(source);

    interpolated_ = true;
}

template<typename Points, typename Values, typename Kernel>
void rbfInterpolationEigen::interpolate 
(
    Points const& points, 
    Values const& values, 
    Kernel kernel, 
    supportType support 
) 
{
    kernel_ = kernel;
    interpolate(points, values, support); 
}

template<typename Point, typename Points, typename Values> 
double rbfInterpolationEigen::evaluate(
    Point const& point, 
    Points const& points, 
    Values const& values, 
    supportType support 
) 
{
    double result = 0.; 

    const auto Npoints = points.size(); 

    if (! interpolated_)
        interpolate(points, values, support);

    for (sizeType i = 0; i < Npoints; ++i)
        result += rbfCoeffs_[i] * 
                  kernel_(Numeric::distance(point, points[i]), supportRadius_); 

    result += rbfCoeffs_[Npoints]; 

    for (char i = 1; i < Npoly_; ++i)
        result += rbfCoeffs_[i + Npoints] * point[i - 1]; 

    return result;
}


} // End namespace RBF

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

