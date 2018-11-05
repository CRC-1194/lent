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

#ifndef rbfInterpolationEigen_H
#define rbfInterpolationEigen_H

#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace RBF 
{

/*---------------------------------------------------------------------------*\
                         Class rbfInterpolationEigen Declaration
\*---------------------------------------------------------------------------*/

class rbfInterpolationEigen final 
{

public:

    // RBF support radii. 
    enum class supportType 
    {
        GLOBAL_SUPPORT, // N-N point connectivity. 
        MAX_CENTROID_RADIUS, 
        MEAN_CENTROID_RADIUS,
        RMS_RADIUS
    }; 

    // Eigen-based type aliases. 
    using realVector = Eigen::VectorXd;  
    using rbfPoint = std::array<double, 3>; 
    using pointVector = std::vector<rbfPoint>;
    using sizeType = decltype(pointVector().size()); 

    // Constructors

        //- Construct null
        rbfInterpolationEigen();

        //- Construct as copy
        rbfInterpolationEigen(const rbfInterpolationEigen&) = default;

        //- Move construct 
        rbfInterpolationEigen(rbfInterpolationEigen&&) = default; 

        //- Construct from points, and the RBF kernel, factorize only.
        template<typename Points, typename Kernel>
        rbfInterpolationEigen
        (
            Points const& points, 
            Kernel kernel, 
            supportType support = supportType::GLOBAL_SUPPORT
        );

        //- Construct from points, and the RBF kernel, solve for values.
        template<typename Points, typename Values, typename Kernel>
        rbfInterpolationEigen
        (
            Points const& points, 
            Values const& values,
            Kernel kernel, 
            supportType support = supportType::GLOBAL_SUPPORT
        );

    //- Destructor
    ~rbfInterpolationEigen() = default;

    // Member Functions

        template<typename Points>
        Eigen::MatrixXd matrix
        (
            Points const& points, 
            supportType support = supportType::GLOBAL_SUPPORT
        );
    
        // Construct and factorize with the initialized Kernel.
        template<typename Points>
        void factorize 
        (
            Points const& points, 
            supportType support = supportType::GLOBAL_SUPPORT
        ); 

        // Construct and factorize with the a new Kernel.
        template<typename Points, typename Values, typename Kernel>
        void factorize 
        (
            Points const& points, 
            Kernel kernel, 
            supportType support = supportType::GLOBAL_SUPPORT
        ); 
    
        //- Construct and solve using the initialized Kernel.
        template<typename Points, typename Values>
        void interpolate
        (
            Points const& points, 
            Values const& values, 
            supportType support = supportType::GLOBAL_SUPPORT
        ); 

        //- Construct and solve using the initialized Kernel.
        template<typename Points, typename Values, typename Kernel>
        void interpolate 
        (
            Points const& points, 
            Values const& values, 
            Kernel kernel,
            supportType support = supportType::GLOBAL_SUPPORT
        ); 

        template<typename Kernel> 
        void setKernel(Kernel kernel)
        {
            kernel_ = kernel;
            interpolated_ = false;
        }

        //- Evaluate the RBF interpolant at a point. 
        template<typename Point, typename Points, typename Values> 
        inline double evaluate
        (
            Point const& evalPoint, 
            Points const& points, 
            Values const& values, 
            supportType support = supportType::GLOBAL_SUPPORT
        );

    // Member Operators
        rbfInterpolationEigen& operator=(const rbfInterpolationEigen&) = default;

        //- Move assign
        rbfInterpolationEigen& operator=(rbfInterpolationEigen&&) = default;

private:

    // 3D linear polynomial coefficient dimension, dim(a0 + a1x + a2y + a3z) = 4. 
    constexpr static int Npoly_= 4; 

    // RBF Matrix factorization.  
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> matrixFactorization_; 

    // RBF interpolant coefficients. 
    Eigen::VectorXd rbfCoeffs_; 

    // RBF kernel function signature.
    using kernelFunction = std::function<double(double, double)>; 

    // RBF kernel function.
    kernelFunction kernel_; 

    // RBF scaling radius for compact support.
    double supportRadius_; 

    // Internal status flags. 
    bool factorized_; 
    bool interpolated_; 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBF

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rbfInterpolationEigen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
