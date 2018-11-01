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
#include "rbfInterpolation.H"

#include <random>
#include <cmath>
#include <algorithm>

#include <iostream>
#include <fstream>
#include "gtest.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace RBF;

TEST(RBF, RANDOM_POINTS)
{
    // Random point generation initialization 
    // Make these variables global.
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0, 1.0);

    // Extract these into functions.
    auto seedPoints = [&dis,&gen](auto& points)
    {
        for (auto& point : points)
        {
            point[0] = dis(gen);  
            point[1] = dis(gen);  
            point[2] = dis(gen);  
        }
    };

    // Extract into a function.
    auto setValues1x2y2z2 = [](const auto& points, auto& values)
    {
        // Set the explicit values equal to 
        // u(x) = 1 + x^2 + y^2, which is the exact solution to the
        // Poisson equation \nabla^2 u = 6.
        const sizeType Nvals = values.size(); 
        assert(Nvals == points.size()); 
        for (sizeType i = 0; i < points.size(); ++i)
            values[i] = 1 + pow(points[i][0], 3) 
                + pow(points[i][1], 3) + pow(points[i][2], 3);
    };

    // RBF interpolation is used within LENT on small stencils.  
    const std::vector<sizeType> testSizes = {16, 32, 64, 128, 256}; 

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
        seedPoints(nodalPointsMap[Npts]);

        // Initialize and set the random nodal values. 
        nodalValuesMap[Npts] = realVector(Npts); 
        setValues1x2y2z2(nodalPointsMap[Npts], nodalValuesMap[Npts]);

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
            auto rbf = make_rbf(
                nodalPoints, 
                nodalValues, 
                [&rbfName](double r, double rs) 
                { return RBF::rbfKernels.at(rbfName)(r, rs, 1.0); }
            );
            
            realVector nodalRbfValues(nodalPoints.size());
            for(sizeType i = 0; i < nodalPoints.size(); ++i)
                nodalRbfValues[i] = rbf.interpolate(nodalPoints[i]);

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
            seedPoints(testPoints);

            realVector testValuesRBF (Ntest); 
            for (sizeType i = 0; i < Ntest; ++i)
                testValuesRBF[i] = rbf.interpolate(testPoints[i]);
            
            realVector testValuesExact (Ntest); 
            setValues1x2y2z2(testPoints, testValuesExact); 

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
