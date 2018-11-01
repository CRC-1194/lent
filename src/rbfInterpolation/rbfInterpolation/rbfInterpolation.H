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

#ifndef rbfInterpolation_H
#define rbfInterpolation_H

#include <Eigen/Dense>
#include <Eigen/Cholesky>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace RBF 
{

using realVector = Eigen::VectorXd;  
using rbfPoint = std::array<double, 3>; 
using pointVector = std::vector<rbfPoint>;
using sizeType = decltype(pointVector().size()); 

/*---------------------------------------------------------------------------*\
                         Class rbfInterpolation Declaration
\*---------------------------------------------------------------------------*/

// Different scaling types used for the rbfInterpolation.
enum class scalingType 
{
    NONE, 
    MAX_CENTROID_RADIUS, 
    MEAN_CENTROID_RADIUS,
    RMS_RADIUS
}; 

template<class Points, class Values>
class rbfInterpolation final
{
    // Private data
    Points const& pointsRef_;  
    Values const& valuesRef_; 

    using ptSizeType = decltype(pointsRef_.size());
    const ptSizeType Npts_; 
    constexpr static char Nbetas_= 4; 

    Eigen::VectorXd x_, b_; 

    using kernelFunction = std::function<double(double, double)>; 
    kernelFunction kernel_; 

    double rs_; 

    // Private Member Functions
    
        //- Construct and solve the linear system of the RBF interpolation.
        void interpolate(scalingType scaling); 

public:

    // Constructors

        //- Construct null
        rbfInterpolation() = delete;

        //- Construct as copy
        rbfInterpolation(const rbfInterpolation&) = default;

        //- Move construct 
        rbfInterpolation(rbfInterpolation&&) = default; 

        //- Construct from points, values and the RBF kernel.
        template<typename Kernel>
        rbfInterpolation
        (
            Points const& points, 
            Values const& values, 
            Kernel kernel, 
            scalingType scaling = scalingType::NONE
        );


    //- Destructor
    ~rbfInterpolation() = default;

    // Member Functions

        //- Evaluate the interpolant (interpolate) at a point. 
        template<typename Point> 
        inline double interpolate(Point const & evalPoint) const;


    // Member Operators
        rbfInterpolation& operator=(const rbfInterpolation&) = delete;

        //- Move assign
        rbfInterpolation&& operator=(rbfInterpolation&&) = delete;

};

// RBF interpolation factory function template.
template<typename Points, typename Values, typename Kernel> 
rbfInterpolation<Points, Values> 
make_rbf(
    Points const& points, 
    Values const& values, 
    Kernel kernel, 
    scalingType scaling = scalingType::MAX_CENTROID_RADIUS
)
{
    return rbfInterpolation<Points,Values>(points, values, kernel, scaling); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBF

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rbfInterpolation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
