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

#include <random>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
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

}

template<typename Points, typename Values> 
class interpolationRBF 
{
    using lambdasContainer  = std::vector<double>;
    using betasContainer = std::array<double, 4>;
    using kernelFunction = std::function<double(double)>;

    Points const& pointsRef_;  
    Values const& valuesRef_; 

    betasContainer betas_ = {0., 0., 0., 0.}; 
    lambdasContainer lambdas_; 
    kernelFunction kernel_; 

    using sizeType = decltype(lambdas_.size());

    void interpolate()
    {
        using namespace Eigen; 
        using namespace Numeric;

        const auto Nlambdas = lambdas_.size(); 
        const auto Nbetas = betas_.size(); 
        const auto Ntotal = Nlambdas + Nbetas;  

        MatrixXd A (Ntotal, Ntotal);
        VectorXd b (Ntotal);

        // Fill the RBF matrix block and set the source term.
        for (sizeType j = 0; j < Nlambdas; ++j)
        {
            b[j] = valuesRef_[j]; 
            for (sizeType i = 0; i < Nlambdas; ++i)
                A(i,j) = kernel_(distance(pointsRef_[i], pointsRef_[j]));
        }

        // Fill the linear polynomial coefficients.
        auto P = A.topRightCorner(Nlambdas, Nbetas);  
        auto Pt = A.bottomLeftCorner(Nbetas, Nlambdas);
        for(sizeType i = 0; i < Nlambdas; ++i)
        {
            P(i, 0) = 1.; 
            Pt(0, i) = 1.;
            for(sizeType k = 1; k < Nbetas; ++k)
            {
                P(i,k) = pointsRef_[i][k-1]; 
                Pt(k,i) = pointsRef_[i][k-1]; 
            }
        }
        A.bottomRightCorner(Nbetas, Nbetas) = MatrixXd::Zero(Nbetas, Nbetas); 
    }

    public:

        template<typename Kernel>
        interpolationRBF(Points const& points, Values const& values, Kernel kernel)
        :
            pointsRef_(points), 
            valuesRef_(values),
            lambdas_(points.size(), 0.), 
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

            for (decltype(lambdas_.size()) i = 0; i < lambdas_.size(); ++i)
                result += lambdas_[i] * 
                          kernel_(Numeric::distance(evalPoint, pointsRef_[i])); 

            result += betas_[0]; 

            for (char i = 1; i < 4; ++i)
                result += betas_[i] * evalPoint[i]; 

            return result;
        }

};

template<typename Points, typename Values, typename Kernel> 
interpolationRBF<Points, Values> 
make_interpolation_rbf(Points const& points, Values const& values, Kernel kernel)
{
    return interpolationRBF<Points,Values>(points, values, kernel); 
}

TEST(RBF, RANDOM_POINTS_UNIT_DOMAIN)
{
    //extern int mainArgc;
    //extern char** mainArgv;

    //int argc = mainArgc;
    //char** argv = mainArgv;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    using realvec = std::vector<double>;  
    using rbfpoint = std::array<double, 3>; 
    using pointvec = std::vector<rbfpoint>;


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
        assert(points.size() == values.size());
        for (decltype(points.size()) i = 0; i < points.size(); ++i)
            values[i] = 1 + pow(points[i][0], 2) + pow(points[i][1], 2);
    };

    //const unsigned int Np = 5; 
    //pointvec nodalPoints (Np, {0,0,0});
    //seedPoints(nodalPoints);
    // Test the RBF matrix assembly for the cube BCC lattice / stencil. 
    pointvec nodalPoints {{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0},
                          {0.,0.,1.},{1.,0.,1.},{1.,1.,1.},{0.,1.,1.}};
    realvec  nodalValues(nodalPoints.size(), 0); 

    pointvec testPoints (nodalPoints.size(), {0,0,0});
    seedPoints(testPoints);
    realvec  testValues (testPoints.size(), 0); 

    setValues1x2y2(nodalPoints, nodalValues);

    auto rbf = make_interpolation_rbf(nodalPoints, nodalValues, [](double r){ return r; });  
    
    for (decltype(testPoints.size()) i = 0; i < testPoints.size(); ++i)
        testValues[i] = rbf.evaluate(testPoints[i]);

    // Compute and report the L_inf, L_2 and L_1 error. 
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
