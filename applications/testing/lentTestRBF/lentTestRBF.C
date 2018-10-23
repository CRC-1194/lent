/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Test application for the interpolation methods used by the LENT method.

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "messageStream.H"

#include <vector>
#include <string>
#include <array>
#include <functional>
#include <algorithm>
#include <random>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

#include <iostream>
#include "gtest.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace std; 

namespace Numeric 
{
    template<typename PointA, typename PointB> 
    double distance(PointA const& pA, PointB const & pB)
    {
        return std::sqrt(pow(pB[0] - pA[0], 2) + 
                         pow(pB[1] - pA[1], 2) + 
                         pow(pB[2] - pA[2], 2));
    }
}

namespace RBF
{

    using rbfFunctionType = std::function<double(double, double, double)>;
    using rbfTableType = map<string, rbfFunctionType>;
    
    // TODO: Some RBF kernels have free parameters, introduce them if necessary.
    rbfTableType rbfKernels = 
    {
        {
            "WENDLAND_C0",
            [](double r, double rs=1.0, double n = 1.) 
            { 
                return std::pow(1. - (r / rs),2); 
            } 
        },
        {
            "WENDLAND_C2",
            [](double r, double rs=1.0, double n = 1.) 
            { 
                const auto eta = r / rs;
                return std::pow(1. - eta, 4) * (4.*eta + 1.); 
            }
        },
        {
            "WENDLAND_C4", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs; 
                return std::pow(1. - eta, 6) * (35.*std::pow(eta, 2)  + 18.*eta + 3.);
            }
        },
        {
            "WENDLAND_C6",
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs;
                return pow((1. - eta), 8) * (32.*std::pow(eta, 3) + 25.*std::pow(eta, 2) + 8.*(eta) +  1.);
            }

        },
        {
            "LINEAR",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return r / rs;
            }
        },
        {
            "CUBIC_SPLINE",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::pow(r / rs, 3);
            }
        },
        {
            "LINEAR_CUBIC_SPLINE",
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs;
                return eta - std::pow(eta, 3);
            }
        },
        {
            "TRUNCATED_POWER",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::pow(1. - (r / rs), n);
            }
        },
        {
            "MULTIQUADRIC",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::sqrt(1. + std::pow(r / rs, 2));
            }
        },
        {
            "INVERSE_MULTIQUADRIC", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                return std::pow(1 + (r / rs), -n);
            }
        },
        {
            "INVERSE_QUADRATIC",
            [](double r, double rs=1.0, double n = 1.) 
            {
                return 1.0  / (1.0 + std::pow(r / rs, 2));
            }
        }, 
        {
            "THIN_PLATE_SPLINE", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs; 
                return eta * eta * std::log(eta + std::numeric_limits<double>::epsilon());
            }
        },
        {
            "GAUSSIAN", 
            [](double r, double rs=1.0, double n = 1.) 
            {
                const auto eta = r / rs; 
                return std::exp(-n * n * eta * eta);
            }
        }
    };


template<typename Points, typename Values> 
class interpolationRBF 
{

    Points const& pointsRef_;  
    Values const& valuesRef_; 

    using ptSizeType = decltype(pointsRef_.size());
    const ptSizeType Npts_; 
    constexpr static char Nbetas_= 4; 

    Eigen::VectorXd x_, b_; 
    Eigen::LDLT<Eigen::MatrixXd> LDLT_; 

    using kernelFunction = std::function<double(double)>; 
    kernelFunction kernel_; 


    void interpolate()
    {
        using namespace Eigen; 
        using namespace Numeric;

        const auto Ntotal = Npts_ + Nbetas_;  

        MatrixXd A (Ntotal, Ntotal);

        // Fill the RBF matrix block and set the source term.
        for (ptSizeType j = 0; j < Npts_; ++j)
        {
            b_[j] = valuesRef_[j]; 
            for (ptSizeType i = 0; i < Npts_ ; ++i)
                // FIXME: 
                // Compute and use the point radius for scaling RBFs. 
                // Check how to incorporate default arguments.
                A(i,j) = kernel_(distance(pointsRef_[i], pointsRef_[j]));
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

        // Solve for the RBF coeffs and store them.
        LDLT_ = A.ldlt(); 
        x_ = LDLT_.solve(b_);
    }

    public:

        template<typename Kernel>
        interpolationRBF(Points const& points, Values const& values, Kernel kernel)
        :
            pointsRef_(points), 
            valuesRef_(values),
            Npts_(pointsRef_.size()), 
            b_(Nbetas_ + Npts_), 
            kernel_(kernel)
        {
            interpolate();
        }

        // TODO: 
        //  - Add a copy-constructor.
        //  - Add a move-constructor.
        //  - Delete the assignment operator. 

        template<typename Point> 
        double evaluate(Point const & evalPoint) const
        {
            double result = 0.; 

            // TODO: Reuse point radius. 
            for (ptSizeType i = 0; i < Npts_; ++i)
                result += x_[i] * 
                          kernel_(Numeric::distance(evalPoint, pointsRef_[i])); 

            result += x_[Npts_]; 

            for (char i = 1; i < 4; ++i)
                result += x_[i + Npts_] * evalPoint[i - 1]; 

            return result;
        }

};

template<typename Points, typename Values, typename Kernel> 
interpolationRBF<Points, Values> 
make_interpolation_rbf(Points const& points, Values const& values, Kernel kernel)
{
    return interpolationRBF<Points,Values>(points, values, kernel); 
}

}

using namespace RBF;

// TODO: RBF, BCC_CUBE

// TODO: RBF, BCC_FVM_CUBE



TEST(RBF, RANDOM_POINTS_UNIT_DOMAIN)
{
    //extern int mainArgc;
    //extern char** mainArgv;

    //int argc = mainArgc;
    //char** argv = mainArgv;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //using realvec = std::vector<double>;  
    using realVector = Eigen::VectorXd;  
    using rbfPoint = std::array<double, 3>; 
    using pointVector = std::vector<rbfPoint>;
    using sizeType = decltype(pointVector().size()); 

    // Random point generation initialization 
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0, 1.0);

    auto seedPoints = [&dis,&gen](auto& points)
    {
        for (auto& point : points)
        {
            point[0] = dis(gen);  
            point[1] = dis(gen);  
            point[2] = dis(gen);  
        }
    };

    auto setValues1x2y2 = [](const auto& points, auto& values)
    {
        // Set the explicit values equal to 
        // u(x) = 1 + x^2 + y^2, which is the exact solution to the
        // Poisson equation \nabla^2 u = 4.
        const sizeType Nvals = values.size(); 
        assert(Nvals == points.size()); 
        for (sizeType i = 0; i < points.size(); ++i)
            values[i] = 1 + pow(points[i][0], 2) + pow(points[i][1], 2);
    };


    // SEEDED POINTS
    // const sizeType Np = 5; 
    // pointVector nodalPoints (Np, {0,0,0});
    // seedPoints(nodalPoints);

    // BCC CUBE 
    // Test the RBF matrix assembly for the cube BCC lattice / stencil. 
    pointVector nodalPoints {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0},
                             {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};

    realVector nodalValues(nodalPoints.size()); 
    setValues1x2y2(nodalPoints, nodalValues);

    // Compute the RBF interpolant.
    auto rbf = make_interpolation_rbf(
        nodalPoints, 
        nodalValues, 
        [](double r) { return RBF::rbfKernels["GAUSSIAN"](r, 1.0, 1.0); }
    );
    
    // Relative errors in the nodal points should be lower than machne epsilon.
    realVector nodalRbfValues(nodalPoints.size());
    for(sizeType i = 0; i < nodalPoints.size(); ++i)
        nodalRbfValues[i] = rbf.evaluate(nodalPoints[i]);

    const auto nodalError = (nodalRbfValues - nodalValues).cwiseAbs();

    auto calcRelError = [](const auto& inputError, const auto& values) 
    {
        realVector relError (inputError);  
        for (decltype(inputError.size()) i = 0; i < inputError.size(); ++i)
            if (relError[i] > 0)
                relError[i] = relError[i] / values[i];
        return relError;
    };

    const auto nodalRelError = calcRelError(nodalError, nodalValues);  
    const auto nodalRelErrorMax = 
        nodalRelError.lpNorm<Eigen::Infinity>();

    const auto epsilon = std::numeric_limits<double>::epsilon();
    EXPECT_TRUE(nodalRelErrorMax <= epsilon)
        << "L1 nodal error = " << nodalRelErrorMax << std::endl
        << "epsilon = " << epsilon; 

    sizeType Ntest = 1000; 
    pointVector testPoints (Ntest, {0,0,0});
    seedPoints(testPoints);

    realVector testValuesRBF (Ntest); 
    for (sizeType i = 0; i < Ntest; ++i)
        testValuesRBF[i] = rbf.evaluate(testPoints[i]);
    
    realVector testValuesExact (Ntest); 
    setValues1x2y2(testPoints, testValuesExact); 

    const auto testError = (testValuesRBF - testValuesExact).cwiseAbs();
    const auto errorL1 = testError.lpNorm<1>() / testError.size(); 
    const auto errorL2 = testError.lpNorm<2>() / testError.size(); 
    const auto errorLinf = testError.lpNorm<Eigen::Infinity>(); 
    const auto errorMean = testError.mean(); 

    //std::cout << errorLinf << "," << errorMean << "," 
        //<< errorL2 << "," << errorL1 << std::endl;
}

int mainArgc;
char** mainArgv;

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    mainArgc = argc;
    mainArgv = argv;
    

    return RUN_ALL_TESTS();

    return 0;
}

// ************************************************************************* //
