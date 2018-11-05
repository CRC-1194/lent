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

Author
    Tomislav Maric maric@mma.tu-darmstadt.de

Description
    Test application for the RBF interpolation over random points. 

    Generates an increased number of random points and interpolates 

    f(x,y,z) = 1 + x^2 + y^2 + z^2 using all functions from rbFunctions.H 

    and checks that the solution of the linear system is under a prescribed
    tolerance at nodal points. 

    TODO: 

    1. Add a map of interesting functions to loop over: 
        1.1 A signed distance to a sphere. 
        1.2 Harmonic potential for the velocity fields used for the advection. 

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "messageStream.H"
#include "rbFunctions.H"
#include "rbfInterpolationEigen.H"

#include <random>
#include <cmath>
#include <algorithm>

#include <iostream>
#include <fstream>
#include "gtest.h"

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

// Harmonic scalar potential in 3D that results in divergence-free
// velocities.
//template<typename Points, typename Values>
//void setValuesShear3D(auto& values, const auto& points, double omega)
//{
    //const sizeType Nvals = values.size(); 
    //assert(Nvals == points.size()); 
    ////for (sizeType i = 0; i < points.size(); ++i)
        ////values[i] = cos(omega * points[i][0]) 
            ///[>cos(points[i][1], n) + pow(points[i][2], n);
//}

//template<int N>
//constexpr decltype(auto) make_set_valuesXnYnZn()
//{
    //return [](auto& values, const auto& points) { return setValues1xnynzn(values, points, N); };
//}

//auto setValuesX2Y2Z2 = make_set_valuesXnYnZn<2>(); 
//auto setValues1x2y2z2 = make_set_values<3>();  
//auto setValues1x2y2z2 = make_set_values<4>();  

using namespace RBF;

TEST(RBF_EIGEN, CLASS_INTERFACE)
{
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0, 1.0);

    using pointVector = rbfInterpolationEigen::pointVector; 
    using realVector = rbfInterpolationEigen::realVector; 

    // TODO: Loop 100 times over the test. 
    const unsigned int testSize = 64; 

    pointVector nodalPoints(testSize); 
    seedPoints(nodalPoints, dis, gen);

    // TODO: Loop over exact value functions. 
    realVector nodalValues(nodalPoints.size()); 
    setValuesXnYnZn(nodalValues, nodalPoints, 2);

    unsigned long int Npoints = nodalPoints.size(); 
    unsigned long int Nvalues = nodalValues.size(); 

    ASSERT_EQ(Npoints, Nvalues);

    // For all RBF kernels.
    for (const auto & nameKernelPair : RBF::rbfKernels)
    {
        const auto rbfName = nameKernelPair.first; 

        auto rbfKernel = [rbfName](double r, double rs) 
        { return RBF::rbfKernels.at(rbfName)(r, rs, 1.0); };

        // Emtpy construction test.
        rbfInterpolationEigen rbfEmpty; 
        // Factorization and sytem solution.
        rbfEmpty.interpolate(nodalPoints, nodalValues, rbfKernel);

        // Construct with factorization. 
        rbfInterpolationEigen rbfFactorized(nodalPoints, rbfKernel);  
        // Solve with given values. 
        // This assumes nodalPoints are those given to the constructor! 
        rbfFactorized.interpolate(nodalPoints, nodalValues);

        // Do all at once: assemble and solve the system.
        rbfInterpolationEigen rbfFull(nodalPoints, nodalValues, rbfKernel);  

        realVector rbfEmptyVals (nodalValues.size()); 
        realVector rbfFactorizedVals (nodalValues.size()); 
        realVector rbfFullVals (nodalValues.size()); 

        for(decltype(nodalValues.size()) i = 0; i < nodalValues.size(); ++i)
        {
            rbfEmptyVals[i] = 
                rbfEmpty.evaluate(nodalPoints[i], nodalPoints, nodalValues);  
            rbfFactorizedVals[i] = 
                rbfFactorized.evaluate(nodalPoints[i], nodalPoints, nodalValues);  
            rbfFullVals[i] = 
                rbfFull.evaluate(nodalPoints[i], nodalPoints, nodalValues);  
        }

        // Solutions constructed by different member functions are equally accurate. 
        const auto eEmpty = rbfEmptyVals - nodalValues;
        const auto lInfEmpty = eEmpty.lpNorm<Eigen::Infinity>();
        EXPECT_LE(lInfEmpty, 1e-10) // Tolerance depends on the linear solver. 
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";
        
        const auto eFactorized = rbfFactorizedVals - nodalValues;
        const auto lInfFactorized = eFactorized.lpNorm<Eigen::Infinity>();
        EXPECT_LE(lInfFactorized, 1e-10) // Tolerance depends on the linear solver. 
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";

        const auto eFull = rbfFullVals - nodalValues;
        const auto lInfFull = eFull.lpNorm<Eigen::Infinity>();
        EXPECT_LE(lInfFull, 1e-10) // Tolerance depends on the linear solver.
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";

        // Solutions constructed by different member functions are the same. 
        const auto eEmptyFactorized = rbfEmptyVals - rbfFactorizedVals;
        const auto lInfEmptyFactorized = eEmptyFactorized.lpNorm<Eigen::Infinity>();
        EXPECT_LE(lInfEmptyFactorized, 1e-15)
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";

        const auto eEmptyFull = rbfEmptyVals - rbfFullVals;
        const auto lInfEmptyFull = eEmptyFull.lpNorm<Eigen::Infinity>();
        EXPECT_LE(lInfEmptyFull, 1e-15)
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";

        const auto eFactorizedFull = rbfFactorizedVals - rbfFullVals;
        const auto lInfFactorizedFull = eFactorizedFull.lpNorm<Eigen::Infinity>();
        EXPECT_LE(lInfFactorizedFull, 1e-15)
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";


        // Interpolation using a copy of the values must produce the same result. 
        auto nodalValuesCopy (nodalValues); 

        rbfEmpty.interpolate(nodalPoints, nodalValuesCopy); 
        rbfFactorized.interpolate(nodalPoints, nodalValuesCopy); 
        rbfFull.interpolate(nodalPoints, nodalValuesCopy); 

        realVector rbfEmptyValsCopy (nodalValues.size()); 
        realVector rbfFactorizedValsCopy (nodalValues.size()); 
        realVector rbfFullValsCopy (nodalValues.size()); 

        for(decltype(nodalValuesCopy.size()) i = 0; i < nodalValuesCopy.size(); ++i)
        {
            rbfEmptyValsCopy[i] = 
                rbfEmpty.evaluate(nodalPoints[i], nodalPoints, nodalValues);  
            rbfFactorizedValsCopy[i] = 
                rbfFactorized.evaluate(nodalPoints[i], nodalPoints, nodalValues);  
            rbfFullValsCopy[i] = 
                rbfFull.evaluate(nodalPoints[i], nodalPoints, nodalValues);  
        }

        const auto eEmptyCopy = rbfEmptyValsCopy - rbfEmptyVals; 
        const auto eEmptyLinf = eEmptyCopy.lpNorm<Eigen::Infinity>(); 
        EXPECT_LE(eEmptyLinf, 1e-15)
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";

        const auto eFactorizedCopy = rbfFactorizedValsCopy - rbfFactorizedVals; 
        const auto eFactorizedLinf = eFactorizedCopy.lpNorm<Eigen::Infinity>(); 
        EXPECT_LE(eFactorizedLinf, 1e-15)
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";

        const auto eFullCopy = rbfFullValsCopy - rbfFullVals; 
        const auto eFullLinf = eFullCopy.lpNorm<Eigen::Infinity>(); 
        EXPECT_LE(eFullLinf, 1e-15)
            << "RBF = " << rbfName << "\n"
            << "Npoints = " << Npoints << "\n";
    }
}

TEST(RBF_EIGEN, RANDOM_POINTS)
{
    // Random point generation initialization 
    // Make these variables global.
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0, 1.0);

    using sizeType = rbfInterpolationEigen::sizeType; 
    using pointVector = rbfInterpolationEigen::pointVector; 
    using realVector = rbfInterpolationEigen::realVector; 


    // RBF interpolation is used within LENT on small stencils.  
    const std::vector<sizeType> testSizes = {16, 32, 64, 128}; 

    // Initialize the nodal point and value data. 
    std::map<int, pointVector> nodalPointsMap; 
    std::map<int, realVector> nodalValuesMap;
    std::map<int, double> hMap;  

    // Compute: random points, exact values and fill distances.
    std::ofstream nodalPointFile("RBF_RANDOM_POINTS_NODAL_VALUES.csv"); 
    nodalPointFile << "N_POINTS,POINT_ID,POINT_X,POINT_Y,POINT_Z, VALUE" << std::endl;
    for (auto Npts : testSizes)
    {
        // Initialize and set the random nodal points. 
        nodalPointsMap[Npts] = pointVector(Npts, {0,0,0}); 
        seedPoints(nodalPointsMap[Npts], dis, gen);

        // Initialize and set the random nodal values. 
        nodalValuesMap[Npts] = realVector(Npts); 
        setValuesXnYnZn(nodalValuesMap[Npts], nodalPointsMap.at(Npts), 2.0);

        // Write nodal points to a file.
        const auto& nodalPoints = nodalPointsMap.at(Npts); 
        const auto& nodalValues = nodalValuesMap.at(Npts);
        decltype(nodalPoints.size()) i = 0; 
        for (const auto& point : nodalPoints)
        {
            nodalPointFile << Npts << "," << i << ","
                << point[0] << "," << point[1] << "," << point[2] 
                << "," << nodalValues[i] << std::endl; 
            ++i; 
        }

        // Compute the fill distance of the point set. 
        hMap[Npts] = Numeric::fillDistance(nodalPoints);
    }
    nodalPointFile.close();

    // Write: test sizes, fill distances, RBF names.
    std::ofstream testSizeFile("RBF_RANDOM_POINTS_TEST_SIZES.csv"); 
    testSizeFile << "TEST_SIZE,FILL_DISTANCE" << std::endl;
    for (const auto testSize : testSizes)
        testSizeFile << testSize << "," << hMap.at(testSize) << std::endl; 
    testSizeFile.close();
    std::ofstream rbfNameFile("RBF_RANDOM_POINTS_RBF_NAMES.csv"); 
    rbfNameFile << "RBF_NAME" << std::endl;
    for (const auto & nameKernelPair : RBF::rbfKernels)
        rbfNameFile << nameKernelPair.first << std::endl;
    rbfNameFile.close();


    std::ofstream errorFile("RBF_RANDOM_POINTS_ERRORS.csv"); 
    errorFile << "RBF,N_POINTS,ERROR_INF,ERROR_RMSE" << std::endl;
    for (const auto & nameKernelPair : RBF::rbfKernels)
    {
        const auto& rbfName = nameKernelPair.first; 
        for (auto Npts : testSizes)
        {
            const auto& nodalPoints = nodalPointsMap.at(Npts); 
            const auto& nodalValues = nodalValuesMap.at(Npts);

            // Compute the RBF interpolant.
            auto rbfKernel = [&rbfName](double r, double rs) 
            { return RBF::rbfKernels.at(rbfName)(r, rs, 1.0); };
            auto rbf = rbfInterpolationEigen(
                nodalPoints, 
                nodalValues, 
                rbfKernel
            );
            
            realVector nodalRbfValues(nodalPoints.size());
            for(sizeType i = 0; i < nodalPoints.size(); ++i)
                nodalRbfValues[i] = rbf.evaluate(nodalPoints[i], nodalPoints, nodalValues);

            const auto nodalErrors = (nodalRbfValues - nodalValues);
            // Use long double for the result to counter error cancellation.
            const long double nodalErrorLinf = nodalErrors.lpNorm<Eigen::Infinity>();
            // TODO: Check this tolerance, all RBFs behave well, except for a few of 
            // them that have larger errors with a larger number of points.
            EXPECT_LE(nodalErrorLinf, 1e-07)
                << "RBF = " << rbfName << std::endl 
                << "Npts = " << Npts << std::endl;

            sizeType Ntest = 10000; // Number of sampling points. 
            pointVector testPoints (Ntest, {0,0,0});
            seedPoints(testPoints, dis, gen);

            realVector testValuesRBF (Ntest); 
            for (sizeType i = 0; i < Ntest; ++i)
                testValuesRBF[i] = rbf.evaluate(testPoints[i], nodalPoints, nodalValues);
            
            realVector testValuesExact (Ntest); 
            setValuesXnYnZn(testValuesExact, testPoints, 2.0); 

            // Simple output: multidimensional indexed CSV for use with pandas.DataFrame.
            const auto testError = (testValuesRBF - testValuesExact);
            const auto Einf =  testError.lpNorm<Eigen::Infinity>(); 
            const auto Ermse = testError.lpNorm<2>() / sqrt(testError.size()); 

            errorFile << rbfName << "," << Npts << "," << Einf << "," 
                << Ermse << std::endl; 
        } 
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
