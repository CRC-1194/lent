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
#include <sstream>
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

// Points of the BCCC RBF cubic stencil
static const pointVector bcccPoints = {
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

// Value Linf is near machine tolerance at nodal points and test points. 
// Gradient Linf is near machine tolerance at nodal points and test points.
template
<
    typename RbfTuple, 
    typename Points, 
    typename RealVector, 
    typename Surface, 
    int N = 0
> 
void test_rbf_plane(
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

    // Are RBF coefficients in the floating point range?
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
            << "Surface = " << surf << std::endl
            << "Point = " << testPoints[pointI] << std::endl;

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
        test_rbf_plane<RbfTuple, Points, RealVector, Surface, N + 1>(
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

        test_rbf_plane<rbfTuple>(
            bccPoints, 
            bccPlaneValues, 
            testPoints, 
            testPlane
        ); 
        
        // BCCC plane test
        realVector bcccPlaneValues (bcccPoints.size()); 
        set_surface_values(bcccPlaneValues, bcccPoints, testPlane);

        test_rbf_plane<rbfTuple>(
            bcccPoints, 
            bcccPlaneValues, 
            testPoints, 
            testPlane
        ); 
    }
}

template
<
    typename RbfTuple, 
    typename Points, 
    typename RealVector, 
    typename Surface, 
    int N = 0
> 
void test_rbf_surface(
    Points const& nodalPoints, 
    RealVector const& nodalValues,
    Surface const& surf, 
    std::string testName, 
    int nTestPoints=3
)
{
    using rbfKernelType = std::tuple_element_t<N, RbfTuple>;
    using rbfInterpolation = rbfInterpolationEigen<rbfKernelType>; 

    // Construct the RBF interpolants, factorize the systems and solve them. 
    rbfInterpolation rbf(nodalPoints, nodalValues);  

    const auto& rbfCoeffs = rbf.coeffs(); 

    // Are RBF coefficients in the floating point range?
    for (decltype(rbfCoeffs.size()) cI = 0; cI < rbfCoeffs.size(); ++cI)
        ASSERT_TRUE((std::numeric_limits<double>::lowest() < rbfCoeffs[cI]) &&
                    (rbfCoeffs[cI] < std::numeric_limits<double>::max()))
            << "RBF coeffs = " << rbfCoeffs << std::endl;

    const double EPS = std::numeric_limits<double>::epsilon();
    
    // Test values at nodal points.
    for (decltype(nodalPoints.size()) pointI = 0; pointI < nodalPoints.size(); ++pointI)
    {
        const auto rbfNodalValue = rbf.value(nodalPoints[pointI], nodalPoints);
        const auto surfValue = surf.value(nodalPoints[pointI]);
        const auto nodalValError = std::abs(rbfNodalValue - surfValue); 
                
        EXPECT_LE(nodalValError, 128 * EPS)
            << rbfKernelType::name() << std::endl
            << "Surface = " << surf << std::endl;
    }

    // Visualize test points, their values, gradients and errors in VTK. 
    // TODO: Extract the legacy VTK output into function.
    const Eigen::Vector3d p0 (0,0,0);
    const double L = 1.0;
    std::ofstream vtks(rbfKernelType::name() + "_" + testName + ".vtk"); 
    vtks << "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET RECTILINEAR_GRID\n"; 
    vtks << "DIMENSIONS " << nTestPoints << " " << nTestPoints << " " << nTestPoints << "\n";
    vtks << "X_COORDINATES " << nTestPoints << " float\n";
    const double xStart = p0[0] + (L / (2*nTestPoints));
    for (int i = 0; i < nTestPoints; ++i)
        vtks << xStart + i*(L / nTestPoints) << " ";
    vtks << "\n";
    vtks << "Y_COORDINATES " << nTestPoints << " float\n";
    const double yStart = p0[1] + (L / (2*nTestPoints));
    for (int i = 0; i < nTestPoints; ++i)
        vtks << yStart + i*(L / nTestPoints) << " ";
    vtks << "\n";
    vtks << "Z_COORDINATES " << nTestPoints << " float\n";
    const double zStart = p0[2] + (L / (2*nTestPoints));
    for (int i = 0; i < nTestPoints; ++i)
        vtks << zStart + i*(L / nTestPoints) << " ";
    vtks << "\n";

    // TODO: LINF value and gradient for every Kernel into CSV, and plot diagrams. 

    vtks << "POINT_DATA " << nTestPoints * nTestPoints * nTestPoints << "\n";
    std::stringstream valss, surfvalss, surfgradss, gradss, gradmagss, valinfss, gradlinfss; 
    valss << "SCALARS rbfval float\nLOOKUP_TABLE DEFAULT\n";
    surfvalss << "SCALARS surfval float\nLOOKUP_TABLE DEFAULT\n";
    gradss << "VECTORS rbfgrad float\n";
    surfgradss << "VECTORS surfgrad float\n";
    gradmagss << "SCALARS rbfgradmag float\nLOOKUP_TABLE DEFAULT\n";
    valinfss << "SCALARS rbfvallinf float\nLOOKUP_TABLE DEFAULT\n";
    gradlinfss << "SCALARS rbfgradlinf float\nLOOKUP_TABLE DEFAULT\n";
    for (int k = 0; k < nTestPoints; ++k)
    {
        for (int j = 0; j < nTestPoints; ++j)
        {
            for (int i = 0; i < nTestPoints; ++i)
            {
                Eigen::Vector3d testPoint(
                    xStart + i * (L / nTestPoints),
                    yStart + j * (L / nTestPoints),
                    zStart + k * (L / nTestPoints) 
                );
                const auto rbfTestValue = rbf.value(testPoint, nodalPoints);
                const auto surfValue = surf.value(testPoint);
                const auto valueError = std::abs(rbfTestValue - surfValue); 

                Eigen::Vector3d rbfTestGrad = rbf.grad(testPoint, nodalPoints); 
                Eigen::Vector3d surfGrad = surf.grad(testPoint); 

                Eigen::Vector3d gradDiff = rbfTestGrad - surfGrad;
                double gradError = std::sqrt(gradDiff.dot(gradDiff));

                valss      << rbfTestValue << " ";
                surfvalss  << surfValue << " ";
                gradss     << rbfTestGrad[0] << " " << rbfTestGrad[1] << " " << rbfTestGrad[2] << " "; 
                surfgradss << surfGrad[0] << " " << surfGrad[1] << " " << surfGrad[2] << " "; 
                gradmagss  << std::sqrt(rbfTestGrad.dot(rbfTestGrad)) << " ";
                valinfss   << valueError << " ";
                gradlinfss << gradError << " ";
            }
                valss      << "\n";
                surfvalss  << "\n";
                gradss     << "\n";
                surfgradss << "\n";
                gradmagss  << "\n";
                valinfss   << "\n";
                gradlinfss << "\n";
        }
    }
    vtks << valss.str() << surfvalss.str() << valinfss.str() 
        << gradss.str() << gradmagss.str() << gradlinfss.str() 
        << surfgradss.str() << "\n"; 

    // Loop over RBF kernels.
    if constexpr (N + 1 < std::tuple_size_v<RbfTuple>) 
        test_rbf_surface<RbfTuple, Points, RealVector, Surface, N + 1>(
                nodalPoints, 
                nodalValues,
                surf,
                testName,
                nTestPoints
            );
}

TEST(RBF_EIGEN, SPHERE_STENCILS)
{
    eigenSphere testSphere(Eigen::Vector3d(0,0,0), 0.750173); 

    // BCC Sphere test
    realVector bccSphereValues(bccPoints.size()); 
    set_surface_values(bccSphereValues, bccPoints, testSphere);

    const int nTestPoints = 50;

    test_rbf_surface<rbfTuple>(
        bccPoints, 
        bccSphereValues, 
        testSphere, 
        "RBF_EIGEN_SPHERE_BCC", 
        nTestPoints 
    ); 

    // BCCC Sphere test
    realVector bcccSphereValues(bcccPoints.size()); 
    set_surface_values(bcccSphereValues, bcccPoints, testSphere);

    test_rbf_surface<rbfTuple>(
        bcccPoints, 
        bcccSphereValues, 
        testSphere,
        "RBF_EIGEN_SPHERE_BCCC", 
        nTestPoints 
    ); 
}

TEST(RBF_EIGEN, ELLIPSOID_STENCILS)
{
    eigenEllipsoid testEllipsoid(1.0 / 3., 0.5, 2./3.); 

    // BCC Ellipsoid test
    realVector bccSphereValues(bccPoints.size()); 
    set_surface_values(bccSphereValues, bccPoints, testEllipsoid);

    const int nTestPoints = 50;

    test_rbf_surface<rbfTuple>(
        bccPoints, 
        bccSphereValues, 
        testEllipsoid, 
        "RBF_EIGEN_ELLIPSOID_BCC", 
        nTestPoints 
    ); 

    // BCCC Ellipsoid test
    realVector bcccSphereValues(bcccPoints.size()); 
    set_surface_values(bcccSphereValues, bcccPoints, testEllipsoid);

    test_rbf_surface<rbfTuple>(
        bcccPoints, 
        bcccSphereValues, 
        testEllipsoid,
        "RBF_EIGEN_ELLIPSOID_BCCC", 
        nTestPoints 
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
