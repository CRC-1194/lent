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

#ifndef rbfInterpolationEigen_H
#define rbfInterpolationEigen_H

#include <Eigen/Dense>
#include "rbfConfig.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBF 
{

/*---------------------------------------------------------------------------*\
                         Class rbfInterpolationEigen Declaration
\*---------------------------------------------------------------------------*/

class rbfInterpolationEigen final 
{

public:

    // Eigen-based type aliases. 
    using realVector = Eigen::VectorXd;  
    using rbfPoint = std::array<double, 3>; 
    using pointVector = std::vector<rbfPoint>;
    using sizeType = decltype(pointVector().size()); 

    using MatrixType = Eigen::Matrix
    <
        double, 
        Eigen::Dynamic, 
        Eigen::Dynamic, 
        Eigen::ColMajor
    >;

    // Constructors

        //- Construct null
        rbfInterpolationEigen();

        //- Construct as copy
        rbfInterpolationEigen(const rbfInterpolationEigen&) = default;

        //- Move construct 
        rbfInterpolationEigen(rbfInterpolationEigen&&) = default; 

        //- Construct the RBF kernel only
        template<typename Kernel>
        rbfInterpolationEigen(
            const Kernel& kernel, 
            const supportType support = supportType::GLOBAL_SUPPORT
        );

        //- Construct from points, and the RBF kernel, factorize only.
        template<typename Points, typename Kernel>
        rbfInterpolationEigen(
            Points const& points, 
            const Kernel& kernel, 
            const supportType support = supportType::GLOBAL_SUPPORT
        );

        //- Construct from points, and the RBF kernel, solve for values.
        template<typename Points, typename Values, typename Kernel>
        rbfInterpolationEigen
        (
            Points const& points, 
            Values const& values,
            const Kernel& kernel,
            const supportType support = supportType::GLOBAL_SUPPORT
        );

    //- Destructor
    ~rbfInterpolationEigen() = default;

    // Member Functions

        // Factorize with the initialized Kernel.
        template<typename Points>
        void factorize(Points const& points); 

        // Factorize with the given Kernel.
        template<typename Points, typename Kernel>
        void factorize(Points const& points, const Kernel& kernel); 

        // Factorize with a new Kernel.
        template<typename Points, typename Values, typename Kernel>
        void factorize (Points const& points, const Kernel&  kernel); 
    
        //- Factorize and solve using the initialized Kernel.
        template<typename Points, typename Values>
        void interpolate (Points const& points, Values const& values); 

        //- Factorize and solve using a new Kernel.
        template<typename Points, typename Values, typename Kernel>
        void interpolate
        (
            Points const& points, 
            Values const& values,
            const Kernel& kernel
        ); 

        //- Evaluate the RBF interpolant at a point. 
        //- Assumes that the RBF has been factorized and interpolated!
        template<typename Point, typename Points, typename Values> 
        inline double evaluate
        (
            Point const& evalPoint, 
            Points const& points, 
            Values const& values
        );

    // Member Operators
        rbfInterpolationEigen& operator=(const rbfInterpolationEigen&) = default;

        //- Move assign
        rbfInterpolationEigen& operator=(rbfInterpolationEigen&&) = default;

private:

    // 3D linear polynomial coefficient dimension, dim(a0 + a1x + a2y + a3z) = 4. 
    constexpr static int Npoly_= 4; 

    // Once it has been factorized, copying makes no sense.

    MatrixType A_; 
    //Eigen::LLT<MatrixType> matrixFactorization_; 
    //Eigen::LDLT<MatrixType> matrixFactorization_; 
    //Eigen::HouseholderQR<MatrixType> matrixFactorization_; 
    Eigen::PartialPivLU<MatrixType> matrixFactorization_; 
    //Eigen::ColPivHouseholderQR<MatrixType> matrixFactorization_; 

    // RBF interpolant coefficients. 
    Eigen::VectorXd rbfCoeffs_; 

    // RBF kernel function signature.
    using kernelFunction = std::function<double(double, double)>; 

    // RBF kernel function.
    kernelFunction kernel_; 

    // RBF support type.
    supportType support_;

    // RBF scaling radius for compact support.
    double supportRadius_; 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RBF
} // End namespace Foam 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rbfInterpolationEigen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
