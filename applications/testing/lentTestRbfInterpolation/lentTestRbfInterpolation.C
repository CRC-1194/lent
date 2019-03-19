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

using namespace Foam;
using namespace RBF;

#include "lentTestRbfInterpolationHelpers.H"

using realVector = Eigen::VectorXd; 
using point = Eigen::Vector3d;  
using pointVector = std::vector<point>; 

// Random number generators.
static std::random_device rd;  
static std::mt19937 gen(rd()); 
static std::uniform_real_distribution<> dis(0, 1.0);
static std::uniform_int_distribution<> dissign(0, 1);

static const std::array<int, 2> signs = {-1, +1}; 

// Points of the BCC RBF cubic stencil
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

// Points of the BCC_FVM RBF cubic stencil
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Tests the emtpy and component constructors and solution member functions
// of the rbfInterpolationEigen<Kernel> class template.
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

// Tests the Linf value errors for nodalValues at nodalPoints and testValues
// at testPoints (tolerance-based).
template
<
    typename RbfTuple, 
    typename Points, 
    typename RealVector, 
    typename Surface, 
    int N = 0
> 
void test_rbf_values_and_grads(
    Points const& nodalPoints, 
    RealVector const& nodalValues,
    Points const& testPoints, 
    Surface const& surf 
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
    
    // Test value and gradient evaluation at nodal points.
    for (decltype(nodalPoints.size()) pointI = 0; pointI < nodalPoints.size(); ++pointI)
    {
        const auto rbfNodalValue = rbf.value(nodalPoints[pointI], nodalPoints);
        const auto surfValue = surf.value(nodalPoints[pointI]);
        const auto nodalError = std::abs(rbfNodalValue - surfValue); 

        ASSERT_LE(nodalError, 8 * EPS)
            << rbfKernelType::name() << std::endl
            << "Surface = " << surf << std::endl;

        Eigen::Vector3d rbfNodalGrad = 
            rbf.template grad<Eigen::Vector3d>(nodalPoints[pointI], nodalPoints); 
        Eigen::Vector3d surfGrad = 
            surf.template grad<Eigen::Vector3d>(nodalPoints[pointI]); 
        Eigen::Vector3d rbfNodalErr = rbfNodalGrad - surfGrad;
        const auto rbfNodalGradErrMag = std::sqrt(rbfNodalErr.dot(rbfNodalErr));
        ASSERT_LE(rbfNodalGradErrMag, 128 * EPS) // Evaluation of RBF gradients near nodal points is unstable.
            << rbfKernelType::name() << std::endl
            << "Surface = " << surf << std::endl
            << "Rbf grad = " << rbfNodalGrad << std::endl
            << "Surf grad = " << surfGrad << std::endl
            << "Error = " << rbfNodalGradErrMag << std::endl
            << "Point = " << nodalPoints[pointI] << std::endl;
    }

    // Gradients and values at test points. 
    for (decltype(testPoints.size()) pointI = 0; pointI < testPoints.size(); ++pointI)
    {
        const auto rbfTestValue = rbf.value(testPoints[pointI], nodalPoints);
        const auto surfValue = surf.value(testPoints[pointI]);
        const auto testError = std::abs(rbfTestValue - surfValue); 

        ASSERT_LE(testError, 8 * EPS)
            << rbfKernelType::name() << std::endl
            << "Surface = " << surf << std::endl;

        Eigen::Vector3d rbfTestGrad = rbf.grad(testPoints[pointI], nodalPoints); 
        Eigen::Vector3d surfGrad = 
            surf.template grad<Eigen::Vector3d>(testPoints[pointI]); 
        Eigen::Vector3d rbfTestErr = rbfTestGrad - surfGrad;
        const auto rbfTestGradErrMag = std::sqrt(rbfTestErr.dot(rbfTestErr));

        ASSERT_LE(rbfTestGradErrMag, 64 * EPS) 
            << "Surface = " << surf << std::endl
            << "Rbf grad = " << rbfTestGrad << std::endl
            << "Surf grad = " << surfGrad << std::endl
            << "Error = " << rbfTestGradErrMag << std::endl
            << "Point = " << testPoints[pointI] << std::endl;
    }
 
    // Loop over RBF kernels.
    if constexpr (N + 1 < std::tuple_size_v<RbfTuple>) 
        test_rbf_values_and_grads<RbfTuple, Points, RealVector, Surface, N + 1>(
                nodalPoints, 
                nodalValues,
                testPoints, 
                surf 
            );
}

TEST(RBF_EIGEN, PLANE_STENCILS)
{
    // Random evaluation points.
    static pointVector testPoints(1000); 
    seedPoints(testPoints, dis, gen); 

    for (int testI = 0; testI < 1000; ++testI)
    {
        eigenPlane testPlane; 
        testPlane.randomize(dis,dissign,gen,signs); 

        // BCC Plane test
        realVector bccPlaneValues (bccPoints.size()); 
        set_surface_values(bccPlaneValues, bccPoints, testPlane);

        test_rbf_values_and_grads<rbfTuple>(
            bccPoints, 
            bccPlaneValues, 
            testPoints, 
            testPlane
        ); 
        
        // BCC_FVM plane test
        realVector bccFvmPlaneValues (bccFvmPoints.size()); 
        set_surface_values(bccFvmPlaneValues, bccFvmPoints, testPlane);

        test_rbf_values_and_grads<rbfTuple>(
            bccFvmPoints, 
            bccFvmPlaneValues, 
            testPoints, 
            testPlane
        ); 
    }
}

TEST(RBF_EIGEN, SPHERE_STENCILS)
{
    // Random evaluation points.
    static pointVector testPoints(1000); 
    seedPoints(testPoints, dis, gen); 

    eigenSphere testSphere(Eigen::Vector3d(0,0,0), 0.50173); 

    // BCC Sphere test
    realVector bccSphereValues(bccPoints.size()); 
    set_surface_values(bccSphereValues, bccPoints, testSphere);

    test_rbf_values_and_grads<rbfTuple>(
        bccPoints, 
        bccSphereValues, 
        testPoints, 
        testSphere 
    ); 

    // BCC_FVM Sphere test
    realVector bccFvmSphereValues(bccFvmPoints.size()); 
    set_surface_values(bccFvmSphereValues, bccFvmPoints, testSphere);

    test_rbf_values_and_grads<rbfTuple>(
        bccFvmPoints, 
        bccFvmSphereValues, 
        testPoints, 
        testSphere 
    ); 
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
