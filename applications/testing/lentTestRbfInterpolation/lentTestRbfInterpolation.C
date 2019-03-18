/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, TU Darmstadt 
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

Description
    A test application for the RBF interpolation. 

Author
    Tomislav Maric maric@mma.tu-darmstadt.de

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "messageStream.H"
#include "rbFunctions.H"
#include "rbfInterpolationEigen.H"

#include <random>
#include <cmath>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <fstream>
#include "gtest.h"

using realVector = Eigen::VectorXd; 
using point = Eigen::Vector3d;  
using pointVector = std::vector<point>; 

static std::random_device rd;  
static std::mt19937 gen(rd()); 
static std::uniform_real_distribution<> dis(0, 1.0);
static std::uniform_int_distribution<> dissign(0, 1);
std::array<int, 2> signs = {-1, +1}; 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Extract these into functions
template<typename Points, typename Distribution, typename Generator>
void seedPoints(
    Points& points, 
    Distribution& dis, 
    Generator& gen
)
{
    for (auto& point : points)
    {
        point[0] = dis(gen);  
        point[1] = dis(gen);  
        point[2] = dis(gen);  
    }
};

template<typename Values, typename Points>
void setValuesXnYnZn(Values& values, Points const& points, double N)
{
    const long unsigned int Npoints = points.size();
    const long unsigned int Nvalues = values.size(); 
    EXPECT_EQ(Npoints, Nvalues);
    for (long unsigned int i = 0; i < Npoints; ++i)
        values[i] = pow(points[i][0], N) 
            + pow(points[i][1], N) + pow(points[i][2], N);
};

using namespace Foam;
using namespace RBF;

template<typename RbfTuple, typename Points, int N = 0> 
void test_rbf_interface(Points const& nodalPoints)
{
    using rbfKernelType = std::tuple_element_t<N, RbfTuple>;
    using rbfInterpolation = rbfInterpolationEigen<rbfKernelType>; 

    realVector nodalValues(nodalPoints.size()); 
    setValuesXnYnZn(nodalValues, nodalPoints, 2);

    unsigned long int Npoints = nodalPoints.size(); 
    unsigned long int Nvalues = nodalValues.size(); 

    ASSERT_EQ(Npoints, Nvalues);

    // TEST empty and full (component constructor)

    // Construct empty. 
    rbfInterpolation rbfEmpty; 

    // Factorize and solve the system.
    rbfEmpty.factorize(nodalPoints); 
    rbfEmpty.solve(nodalPoints, nodalValues);

    // Construct, factorize and solve. 
    rbfInterpolation rbfFull(nodalPoints, nodalValues);  

    realVector rbfEmptyVals (nodalValues.size()); 
    realVector rbfFullVals (nodalValues.size()); 

    for(decltype(nodalValues.size()) i = 0; i < nodalValues.size(); ++i)
    {
        rbfEmptyVals[i] = 
            rbfEmpty.value(nodalPoints[i], nodalPoints);  
        rbfFullVals[i] = 
            rbfFull.value(nodalPoints[i], nodalPoints);  
    }

    const auto eEmpty = rbfEmptyVals - nodalValues;
    const auto lInfEmpty = eEmpty.lpNorm<Eigen::Infinity>();
    EXPECT_LE(lInfEmpty, 1e-10) // Tolerance depends on the linear solver. 
        << "RBF = " << rbfKernelType::name() << "\n"
        << "Npoints = " << Npoints << "\n";

    const auto eFull = rbfFullVals - nodalValues;
    const auto lInfFull = eFull.lpNorm<Eigen::Infinity>();
    EXPECT_LE(lInfFull, 1e-10) // Tolerance depends on the linear solver.
        << "RBF = " << rbfKernelType::name() << "\n"
        << "Npoints = " << Npoints << "\n";

    const auto eEmptyFull = rbfEmptyVals - rbfFullVals;
    const auto lInfEmptyFull = eEmptyFull.lpNorm<Eigen::Infinity>();
    EXPECT_LE(lInfEmptyFull, 1e-15)
        << "RBF = " << rbfKernelType::name() << "\n"
        << "Npoints = " << Npoints << "\n";

    // Interpolation using a copy of the values must produce the same result. 
    auto nodalValuesCopy (nodalValues); 

    rbfEmpty.solve(nodalPoints, nodalValuesCopy); 
    rbfFull.solve(nodalPoints, nodalValuesCopy); 

    realVector rbfEmptyValsCopy (nodalValues.size()); 
    realVector rbfFullValsCopy (nodalValues.size()); 

    for(decltype(nodalValuesCopy.size()) i = 0; i < nodalValuesCopy.size(); ++i)
    {
        rbfEmptyValsCopy[i] = 
            rbfEmpty.value(nodalPoints[i], nodalPoints);  
        rbfFullValsCopy[i] = 
            rbfFull.value(nodalPoints[i], nodalPoints);  
    }

    const auto eEmptyCopy = rbfEmptyValsCopy - rbfEmptyVals; 
    const auto eEmptyLinf = eEmptyCopy.lpNorm<Eigen::Infinity>(); 
    EXPECT_LE(eEmptyLinf, 1e-15)
        << "RBF = " << rbfKernelType::name() << "\n"
        << "Npoints = " << Npoints << "\n";


    const auto eFullCopy = rbfFullValsCopy - rbfFullVals; 
    const auto eFullLinf = eFullCopy.lpNorm<Eigen::Infinity>(); 
    EXPECT_LE(eFullLinf, 1e-15)
        << "RBF = " << rbfKernelType::name() << "\n"
        << "Npoints = " << Npoints << "\n";

    // - Call yourself using the next RBF kernel in the tuple. 
    if constexpr (N + 1 < std::tuple_size_v<RbfTuple>)
        test_rbf_interface<RbfTuple, Points, N+1>(nodalPoints);
}

TEST(RBF_EIGEN, CLASS_INTERFACE)
{
    const unsigned int maxTestSize = 32; 
    pointVector nodalPoints(maxTestSize); 
    seedPoints(nodalPoints, dis, gen);
    test_rbf_interface<rbfTuple>(nodalPoints); 
}

static const pointVector bccPoints = {
        point{0,0,0}, // Cube corner points.
        point{1,0,0},
        point{1,1,0},
        point{0,1,0},
        point{0,0,1},
        point{1,0,1},
        point{1,1,1},
        point{0,1,1},
        point{0.5,0.5,0.5} // Cube centroid.
    };

static const pointVector bccFvmPoints = {
        point{0,0,0}, // Cube corner points.
        point{1,0,0},
        point{1,1,0},
        point{0,1,0},
        point{0,0,1},
        point{1,0,1},
        point{1,1,1},
        point{0,1,1},
        point{0.5,0.5,0.5}, // Cube centroid.
        point{-0.5,0.5,0.5}, // Neighbor cube centroids.
        point{0.5,-0.5,0.5},
        point{1.5,0.5,0.5},
        point{0.5,1.5,0.5},
        point{0.5,0.5,-0.5},
        point{0.5,0.5,1.5},
    };

template<typename RbfTuple, typename Points, typename RealVector, typename Plane, int N = 0> 
void test_rbf_values_and_grads(
    Points const& nodalPoints, 
    RealVector const& nodalValues,
    Points const& testPoints, 
    Plane const& testPlane
)
{
    using rbfKernelType = std::tuple_element_t<N, RbfTuple>;
    using rbfInterpolation = rbfInterpolationEigen<rbfKernelType>; 

    // Construct the RBF interpolants, factorize the systems and solve them. 
    rbfInterpolation rbf(nodalPoints, nodalValues);  

    const auto& rbfCoeffs = rbf.coeffs(); 

    for (decltype(rbfCoeffs.size()) cI = 0; cI < rbfCoeffs.size(); ++cI)
        ASSERT_TRUE((std::numeric_limits<double>::lowest() < rbfCoeffs[cI]) &&
                    (rbfCoeffs[cI] < std::numeric_limits<double>::max()))
            << "RBF coeffs = " << rbfCoeffs << std::endl;

    const double EPS = std::numeric_limits<double>::epsilon();
    for (decltype(nodalPoints.size()) pointI = 0; pointI < nodalPoints.size(); ++pointI)
    {
        const auto rbfNodalValue = rbf.value(nodalPoints[pointI], nodalPoints);
        const auto testPlaneValue = testPlane.value(nodalPoints[pointI]);
        const auto nodalError = std::abs(rbfNodalValue - testPlaneValue); 

        ASSERT_LE(nodalError, 8 * EPS)
            << rbfKernelType::name() << std::endl
            << "Plane point = " << testPlane.point() << std::endl
            << "Plane normal = " << testPlane.normal() << std::endl;

        Eigen::Vector3d rbfNodalGrad = rbf.grad(nodalPoints[pointI], nodalPoints); 
        Eigen::Vector3d rbfNodalErr = rbfNodalGrad - testPlane.normal(); 
        const auto rbfNodalGradErrMag = std::sqrt(rbfNodalErr.dot(rbfNodalErr));

        ASSERT_LE(rbfNodalGradErrMag, 128 * EPS) // Evaluation of RBF gradients near nodal points is unstable.
            << rbfKernelType::name() << std::endl
            << "Plane point = " << testPlane.point() << std::endl
            << "Plane normal = " << testPlane.normal() << std::endl
            << "Rbf grad = " << rbfNodalGrad << std::endl
            << "Error = " << rbfNodalGradErrMag << std::endl;
    }
 
    for (decltype(testPoints.size()) pointI = 0; pointI < testPoints.size(); ++pointI)
    {
        const auto rbfValue = rbf.value(testPoints[pointI], nodalPoints);
        const auto planeValue = testPlane.value(testPoints[pointI]);
        const auto testError = std::abs(rbfValue - planeValue); 

        ASSERT_LE(testError, 8 * EPS)
            << rbfKernelType::name() << std::endl;
    }

    // GENERIC TEST LOOP  
    // - Call yourself for the next RBF kernel or stop.
    if constexpr (N + 1 < std::tuple_size_v<RbfTuple>) 
        test_rbf_values_and_grads<RbfTuple, Points, RealVector, Plane, N + 1>(
                nodalPoints, 
                nodalValues,
                testPoints, 
                testPlane
            );
}

template<typename Point, typename DIS, typename GEN>
Point rand_point(DIS& dis, GEN& gen)
{
    Point pt; 
    pt[0] = dis(gen);
    pt[1] = dis(gen);
    pt[2] = dis(gen);
    return pt;
}

class plane
{
    Eigen::Vector3d point_; 
    Eigen::Vector3d normal_;

    public: 

        plane()
            :
                point_(), 
                normal_()
        {}

        template<typename DIS, typename DISSIGN, typename GEN, typename SIGNS>
        void randomize(DIS& dis, DISSIGN& dissign, GEN& gen, SIGNS const& signs)
        {
            point_ = rand_point<Eigen::Vector3d>(dis, gen);
            normal_ = rand_point<Eigen::Vector3d>(dis, gen);
            normal_[0] = signs[dissign(gen)] * normal_[0];
            normal_[1] = signs[dissign(gen)] * normal_[1];
            normal_[2] = signs[dissign(gen)] * normal_[2];
            normal_.normalize();  
        }

        template<typename Vector>
        decltype(auto) value(Vector&& evalPoint) const
        {
            return (evalPoint - point_).dot(normal_); 
        }

        decltype(auto) normal() const
        {
            return normal_;
        }

        decltype(auto) point() const
        {
            return point_;
        }
};

template<typename Values, typename Points, typename Plane> 
void set_plane_values(
    Values& vals, 
    Points const& points, 
    Plane const& plane
)
{
    for (decltype(points.size()) pI = 0; pI < points.size(); ++pI)
        vals[pI] = plane.value(points[pI]);
}

TEST(RBF_EIGEN, PLANE_SIGNED_DISTANCE)
{

    pointVector testPoints(1000); 
    seedPoints(testPoints, dis, gen); 

    for (int testI = 0; testI < 1000; ++testI)
    {
        plane testPlane; 
        testPlane.randomize(dis,dissign,gen,signs); 

        // BCC Stencil
        
        // BCC Plane test
        realVector bccPlaneValues (bccPoints.size()); 
        set_plane_values(bccPlaneValues, bccPoints, testPlane);

        test_rbf_values_and_grads<rbfTuple>(
            bccPoints, 
            bccPlaneValues, 
            testPoints, 
            testPlane
        ); 

        // BCC_FVM "1 + x**2 + y**2 + z**2" test 
        
        //realVector bcc1x2y2z2Values (bccPoints.size()); 
        //test_rbf_values_and_grads<rbfTuple>(
            //bccPoints, 
            //bcc1x2y2z2Values, 
            //testPoints, 
            //testPlane
        //); 

        // END BCC Stencil

        // BCC_FVM Stencil
        
        // BCC_FVM plane test
        //realVector bccFvmPlaneValues (bccFvmPoints.size()); 
        //set_plane_values(bccFvmPlaneValues, bccFvmPoints, testPlane);

        //test_rbf_values_and_grads<rbfTuple>(
            //bccFvmPoints, 
            //bccFvmPlaneValues, 
            //testPoints, 
            //testPlane
        //); 
        
        // BCC_FVM 1 + x**2 + y**2 + z**2 test 
        //realVector bccFvmPlaneValues (bccFvmPoints.size()); 
        //set_plane_values(bccFvmPlaneValues, testPoints, testPlane);

        //test_rbf_values_and_grads<rbfTuple>(
            //bccFvmPoints, 
            //bccFvmPlaneValues, 
            //testPoints, 
            //testPlane
        //); 
        
        // END BCC_FVM Stencil
    }
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
