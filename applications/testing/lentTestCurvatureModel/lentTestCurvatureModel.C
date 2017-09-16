/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by orignal authors.
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

Application
    lentTestCurvatureModels 

Description
    Test LENT curvature models against the exact curvature. 

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "lentMethod.H"
#include "volMesh.H"
#include "map.H"
#include <fstream>

#include "analyticalSurface.H"
#include "errorMetrics.H"
#include <map>
#include <random>
#include <vector>

using namespace FrontTracking;

// Test if file is empty: for an empty file the current position in a stream in
// append-mode is 0
bool fileIsEmpty(std::fstream& file)
{
    if (file.tellg() == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

scalar frontRadius(const triSurfaceFront& front)
{
    scalar R = 0.0;
    vector centre{4, 4, 4};

    const auto& points = front.points();

    forAll(points, I)
    {
        R += mag(points[I] - centre);
    }

    return R/points.size();
}

std::vector<scalar> signedDistanceDeviation(const volScalarField& signedDistance, const analyticalSurface& surface)
{
    std::vector<scalar> deviations{};

    auto exactDistance = 0.0;
    const auto& C = signedDistance.mesh().C();

    forAll(signedDistance, I)
    {
        if (mag(signedDistance[I]) < GREAT)
        {
            exactDistance = surface.signedDistance(C[I]);
            deviations.push_back(signedDistance[I] - exactDistance);
        }
    }

    return deviations;
}

std::vector<scalar> frontDev(const triSurfaceFront& front, const analyticalSurface& surface)
{
    std::vector<scalar> deviations{};

    const auto& points = front.points();

    deviations.resize(points.size());

    forAll(points, I)
    {
        deviations[I] = surface.signedDistance(points[I]);
    }

    return deviations;
}

void correctAndPerturbDistances(volScalarField& signedDistance, const analyticalSurface& surface, const scalar& noiseLevel)
{
    const auto& C = signedDistance.mesh().C();

    // Add noise to signed distances
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<> dis{-noiseLevel, noiseLevel};

    forAll(C, I)
    {
        if (mag(signedDistance[I]) < GREAT)
        {
            signedDistance[I] = surface.signedDistance(C[I]) + dis(gen);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    std::fstream errorFile; 

    errorFile.open("curvatureErrors.dat", std::ios_base::app);

    triSurfaceFront front(
        IOobject(
            "front",
            "front",
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    lentMethod lent(front, mesh);

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    runTime++;

    lent.reconstructFront(front, signedDistance, pointSignedDistance);

    lent.calcSignedDistances(
        signedDistance,
        pointSignedDistance,
        searchDistanceSqr,
        pointSearchDistanceSqr,
        front
    );

    dimensionedScalar h = max(mag(mesh.delta())); 

    // ------- Begin modifications
    Info << "Radius after reconstruction: " << frontRadius(front) << '\n';

    // Set exact signed distance and optionally perturb it
    auto analyticalSurfaceTmp = tmp<analyticalSurface>{analyticalSurface::New(lent.dict().subDict("analyticalSurface"))};
    const auto& surface = analyticalSurfaceTmp.ref();

    // Record deviations in signed distance
    auto deviations = signedDistanceDeviation(signedDistance, surface);
    errorMetrics metrics{deviations};

    Info << "Distance: mean deviation = " << metrics.arithmeticMeanError()
         << ", median = " << metrics.medianError()
         << ", std deviation = " << metrics.standardDeviation()
         << '\n';

    auto distribution = metrics.errorDistribution(100);
    auto normalizedDistribution = metrics.errorDistributionNormalized(100);

    std::fstream distFile;

    distFile.open("distanceErrorDistribution_" + std::to_string(h.value()) + ".dat", std::ios_base::out);
    distFile << "# dev | nPoints | fraction\n";

    for (const auto& entry : distribution)
    {
        distFile << entry.first << ' ' << entry.second << ' ' << normalizedDistribution[entry.first] << std::endl;
    }

    // Record deviations in the front points after reconstruction
    auto pointDev = frontDev(front, surface);
    errorMetrics frontMetrics{pointDev};

    Info << "Front: mean deviation = " << frontMetrics.arithmeticMeanError()
         << ", median = " << frontMetrics.medianError()
         << ", std deviation = " << frontMetrics.standardDeviation()
         << '\n';

    auto frontDist = frontMetrics.errorDistribution(100);
    auto frontDistNorm = frontMetrics.errorDistributionNormalized(100);

    std::fstream frontFile;
    frontFile.open("frontErrorDistribution_" + std::to_string(h.value()) + ".dat", std::ios_base::out);
    frontFile << "# dev | nPoints | fraction\n";

    for (const auto& entry : frontDist)
    {
        frontFile << entry.first << ' ' << entry.second << ' ' << frontDistNorm[entry.first] << std::endl;
    }

    // Write statistic metrics
    std::fstream statisticsFile;
    statisticsFile.open("distanceStatistics.dat", std::ios_base::app);

    if (fileIsEmpty(statisticsFile))
    {
        statisticsFile << "# mesh spacing | distance error mean | distance error std dev | front error mean | front errot std dev\n";
    }
    statisticsFile << h.value() << ' ' << metrics.arithmeticMeanError() << ' '
                   << metrics.standardDeviation() << ' ' << frontMetrics.arithmeticMeanError() << ' '
                   << frontMetrics.standardDeviation() << '\n';

    //correctAndPerturbDistances(signedDistance, surface, 0.0);
    // -------- End modifications

    const dictionary& lentDict = lent.dict(); 

    tmp<frontCurvatureModel> exactCurvatureModelTmp = frontCurvatureModel::New(lentDict.subDict("exactCurvatureModel")); 
    const frontCurvatureModel& exactCurvatureModel = exactCurvatureModelTmp.ref(); 
    tmp<volScalarField> cellCurvatureExactTmp = exactCurvatureModel.cellCurvature(mesh,front);  
    volScalarField& exactCurvature = cellCurvatureExactTmp.ref();  

    const auto& surfaceTensionDict = lentDict.subDict("surfaceTensionForceModel");
    tmp<frontCurvatureModel> numericalCurvatureModelTmp = frontCurvatureModel::New(surfaceTensionDict.subDict("curvatureModel"));  
    const frontCurvatureModel& numericalCurvatureModel = numericalCurvatureModelTmp.ref(); 
    tmp<volScalarField> numericalCurvatureTmp = numericalCurvatureModel.cellCurvature(mesh,front);  
    volScalarField& numericalCurvature = numericalCurvatureTmp.ref();  
    numericalCurvature.rename("numericalCurvature"); 
    numericalCurvature.writeOpt() = IOobject::AUTO_WRITE; 

    // Use the gradient of the markerField as filter since this gives the cells
    // where the curvature is used when the CSF model is used for surface tension
    dimensionedScalar dSMALL("SMALL", pow(dimLength,-1), SMALL);
    volScalarField onesFilter = pos(mag(fvc::grad(markerField)) - 1e4*dSMALL);
    onesFilter.rename("ones"); 
    onesFilter.writeOpt() = IOobject::AUTO_WRITE; 
    numericalCurvature *= onesFilter; 
    exactCurvature *= onesFilter; 

    Info << "\nMaximum numerical curvature: " << max(numericalCurvature).value()
         << "\nMinimum numerical curvature: " << min(numericalCurvature).value() << endl; 

    volScalarField LinfField ("LinfCurvatureErr", mag(exactCurvature - numericalCurvature)/(dSMALL + mag(exactCurvature))); 
    LinfField.writeOpt() = IOobject::AUTO_WRITE; 
    dimensionedScalar Linf = max(LinfField);

    // Compute and write error distribution of curvature
    std::vector<scalar> curvatureErrors{};

    forAll(onesFilter, I)
    {
        if (onesFilter[I] > 0.0)
        {
            curvatureErrors.push_back((exactCurvature[I] - numericalCurvature[I])/(SMALL + mag(exactCurvature[I])));
        }
    }

    errorMetrics curvatureMetrics{curvatureErrors};

    Info << "Curvature: mean deviation = " << curvatureMetrics.arithmeticMeanError()
         << ", median = " << curvatureMetrics.medianError()
         << ", std deviation = " << curvatureMetrics.standardDeviation()
         << '\n';

    std::fstream curvatureFile;

    curvatureFile.open("curvatureErrorDistribution.dat", std::ios_base::out);
    curvatureFile << "# dev | nCells | fraction\n";
    auto absDistribution = curvatureMetrics.errorDistribution(50);
    auto relDistribution = curvatureMetrics.errorDistributionNormalized(50);

    for (const auto& entry : absDistribution)
    {
        curvatureFile << entry.first << ' ' << entry.second << ' ' << relDistribution[entry.first] << '\n';
    }

    // Write results
    if (fileIsEmpty(errorFile))
    {
        errorFile << "# mesh spacing | Linf relative | exact model | numerical model\n";
    }

    errorFile << h.value() << " " << Linf.value() << " " << exactCurvatureModel.type()
              << " " << numericalCurvatureModel.type() << std::endl;

    front.write();
    runTime.writeNow();

    Info <<"Maximal relative curvature error = " << Linf.value() << '\n' << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
