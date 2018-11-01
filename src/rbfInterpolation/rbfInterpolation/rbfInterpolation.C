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
    Foam::rbfInterpolation

Description
    RBF interpolation extended with a linar polynomial. 

    The interpolant is constructed using a colllection of points and a
    collection of values. 

    The interpolation uses the Eigen library is to find the solution to the 
    linear system. LDLT decomposition is used for the system solution. 

    It is not applicable for a large number of points/values. 

SourceFiles
    rbfInterpolation.H
    rbfInterpolation.C

\*---------------------------------------------------------------------------*/

#include "rbfInterpolation.H"
#include "rbfGeometry.H"

#include <limits>

namespace RBF 
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<typename Points, typename Values> 
template<typename Kernel>
rbfInterpolation<Points, Values>::rbfInterpolation
(
    Points const& points, 
    Values const& values, 
    Kernel kernel, 
    scalingType scaling
)
:
            pointsRef_(points), 
            valuesRef_(values),
            Npts_(pointsRef_.size()), 
            b_(Nbetas_ + Npts_), 
            kernel_(kernel), 
            rs_(1.0)
{
    interpolate(scaling); 
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<typename Points, typename Values> 
void rbfInterpolation<Points, Values>::interpolate(scalingType scaling)
{
        using namespace Eigen; 
        using namespace Numeric;

        const auto Ntotal = Npts_ + Nbetas_;  

        MatrixXd A (Ntotal, Ntotal);

        // Compute the scaling RBF radius as the maximal radius between a point 
        // p and the centroid of the set of points in O(2n). 
        if (scaling == scalingType::MAX_CENTROID_RADIUS)
            rs_ = maxCentroidRadius(pointsRef_);  

        if (scaling == scalingType::MEAN_CENTROID_RADIUS)
            rs_ = meanCentroidRadius(pointsRef_);  

        if (scaling == scalingType::RMS_RADIUS)
            rs_ = rmsRadius(pointsRef_);  

        // Fill the RBF matrix block and set the source term.
        for (ptSizeType j = 0; j < Npts_; ++j)
        {
            b_[j] = valuesRef_[j]; 
            for (ptSizeType i = 0; i < Npts_ ; ++i)
                A(i,j) = kernel_(distance(pointsRef_[i], pointsRef_[j]), rs_);
        }
        b_.tail(Nbetas_) = VectorXd::Zero(Nbetas_); 

        // Fill the linear polynomial coefficients.
        auto P = A.topRightCorner(Npts_, Nbetas_);  
        auto Pt = A.bottomLeftCorner(Nbetas_, Npts_);
        for(ptSizeType i = 0; i < Npts_; ++i)
        {
            P(i, 0) = 1.; 
            Pt(0, i) = 1.;
            for(char k = 1; k < Nbetas_; ++k)
            {
                P(i,k) = pointsRef_[i][k-1]; 
                Pt(k,i) = pointsRef_[i][k-1]; 
            }
        }
        A.bottomRightCorner(Nbetas_, Nbetas_) = MatrixXd::Zero(Nbetas_, Nbetas_); 

        // Solve for the RBF interpolation coeffs.
        x_ = A.colPivHouseholderQr().solve(b_);
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


template<typename Points, typename Values> 
template<typename Point> 
double rbfInterpolation<Points, Values>::interpolate(Point const& evalPoint) const
{
    double result = 0.; 

    for (ptSizeType i = 0; i < Npts_; ++i)
        result += x_[i] * 
                  kernel_(Numeric::distance(evalPoint, pointsRef_[i]), rs_); 

    result += x_[Npts_]; 

    for (char i = 1; i < 4; ++i)
        result += x_[i + Npts_] * evalPoint[i - 1]; 

    return result;
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

} // End namespace RBF

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //

