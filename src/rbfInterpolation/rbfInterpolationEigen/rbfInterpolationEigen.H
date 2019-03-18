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

#ifndef rbfInterpolationEigen_H
#define rbfInterpolationEigen_H

#include <Eigen/Dense>
#include "rbfConfig.H"
#include <memory>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam { namespace RBF { 

/*---------------------------------------------------------------------------*\
                         Class rbfInterpolationEigen Declaration
\*---------------------------------------------------------------------------*/

template<typename Kernel>
class rbfInterpolationEigen final 
{

public:

    // Extract this as a template parameter?
    using MatrixType = Eigen::Matrix
    <
        double, 
        Eigen::Dynamic, 
        Eigen::Dynamic, 
        Eigen::ColMajor, 
        37, // Max row size for a regular dodecahedron cell.
        37  // Max column size for a regular dodecahedron cell.
    >;

    // Constructors

        //- Construct null
        rbfInterpolationEigen();

        //- Construct as copy
        rbfInterpolationEigen(const rbfInterpolationEigen&) = default;

        //- Move construct 
        rbfInterpolationEigen(rbfInterpolationEigen&&) = default; 

        //- Construct 
        rbfInterpolationEigen (const supportType support);

        //- Construct and factorize. 
        template<typename Points>
        rbfInterpolationEigen
        (
            Points const& points, 
            const supportType support = supportType::GLOBAL_SUPPORT
        );

        //- Construct, factorize, solve. 
        template<typename Points, typename Values>
        rbfInterpolationEigen
        (
            Points const& points, 
            Values const& values, 
            const supportType support = supportType::GLOBAL_SUPPORT
        );

    //- Destructor
    ~rbfInterpolationEigen() = default;

    // Member Functions

        // Factorize with the given Kernel.
        template<typename Points>
        void factorize(Points const& points); 
    
        //- Factorize and solve using the initialized Kernel.
        template<typename Points, typename Values>
        void solve(Points const& points, Values const& values); 

        //- Evaluate the RBF interpolant at a point. 
        //- Assumes that the RBF has been factorized and interpolated!
        template<typename Point, typename Points> 
        inline double value(Point const& evalPoint, Points const& points) const;

        //- Evaluate the RBF gradient at a point. 
        //- Assumes that the RBF has been factorized and interpolated!
        template<typename Vector=Eigen::Vector3d, typename Points> 
        inline Vector grad(Vector const& evalPoint, Points const& points) const;

    // Member Operators
        rbfInterpolationEigen& operator=(const rbfInterpolationEigen&) = default;

        //- Move assign
        rbfInterpolationEigen& operator=(rbfInterpolationEigen&&) = default;

        std::string name() const
        {
            return kernel_.name();
        }

        const Eigen::VectorXd& coeffs() const
        {
            return rbfCoeffs_;
        }

        decltype(auto) nCoeffs() const
        {
            return rbfCoeffs_.size();
        }

private:

    // 3D linear polynomial coefficient dimension, dim(a0 + a1x + a2y + a3z) = 4. 
    constexpr static int Npoly_= 4; 

    const Kernel kernel_; 

    const supportType support_; 

    double supportRadius_; 

    MatrixType matrix_; 

    Eigen::PartialPivLU<MatrixType> matrixFactorization_; 

    Eigen::VectorXd source_; 

    Eigen::VectorXd rbfCoeffs_; 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}} // End namespace Foam::RBF

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rbfInterpolationEigen.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
