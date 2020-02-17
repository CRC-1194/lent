/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "lentSubalgorithmTest.H"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

namespace Foam {
namespace FrontTracking {

    using clock = std::chrono::high_resolution_clock;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
std::string lentSubalgorithmTest::assembleFilePath() const
{
    std::string dataFilePath = mesh_.time().rootPath() + "/"
                       + mesh_.time().globalCaseName() + "/";

    return dataFilePath;
}

std::string lentSubalgorithmTest::printFoamVector(const vector& v) const
{
    // Enable streaming of Foam::vector into std::fstream since
    // the stream operator "operator<<(std::fstream&, Foam::vector&)"
    // is not defined
    std::ostringstream streamedVector;

    streamedVector << "(" << v.x() << " " << v.y() << " " << v.z() << ")";

    return std::string{streamedVector.str()};
}

void lentSubalgorithmTest::setNextRun() const
{
    auto& runTime = const_cast<Time&>(mesh_.time());
    runTime++;
}

void lentSubalgorithmTest::writeFields() const
{
    if (writeFields_)
    {
        auto& runTime = const_cast<Time&>(mesh_.time());
        runTime.write();
        frontRef().write();
    }
}

void lentSubalgorithmTest::writeMetrics(const word& fileName, const label runNumber) const
{
    // In the first run create result file and write its header
    if (runNumber == 0)
    {
        std::fstream dataFile(assembleFilePath() + fileName, std::ios_base::out);

        dataFile << metricHeader(scalarMetrics_)
                 << metricHeader(vectorMetrics_) << std::endl;

        dataFile.close();
    }

    std::fstream dataFile(assembleFilePath() + fileName, std::ios_base::app);

    // For now: hardcoded precision using scientific notation
    dataFile << std::scientific << std::setprecision(10);

    for (auto const& metricField: scalarMetrics_)
    {
        dataFile << metricField.second[runNumber] << separator_;
    } 

    for (auto const& metricField: vectorMetrics_)
    {
        dataFile << printFoamVector(metricField.second[runNumber]) << separator_;
    } 
    
    dataFile << std::endl;

    dataFile.close();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void lentSubalgorithmTest::addMeasure(const std::string& metricName, const scalar& value)
{
    scalarMetrics_[metricName].push_back(value);
}

void lentSubalgorithmTest::addMeasure(const std::string& metricName, const vector& value)
{
    vectorMetrics_[metricName].push_back(value);
}

void lentSubalgorithmTest::computeExactSignedDistances() const
{
    const auto& frontSurface = surfaceTmp_.ref();

    // cell-centred distances
    auto& signedDistance = lookupSignedDistance();
    frontSurface.setDistance(signedDistance);
    
    // cell-corner distances
    auto& pointSignedDistance = lookupPointSignedDistance();
    frontSurface.setDistance(pointSignedDistance);
}

void lentSubalgorithmTest::computeExactSignedDistancesNarrowBand() const
{
    // Only compute the exact signed distance in the narrow band
    const auto& frontSurface = surfaceTmp_.ref();

    // cell-centred distances
    auto& signedDistance = lookupSignedDistance();
    const auto& C = mesh().C();

    forAll(signedDistance, I)
    {
        if (mag(signedDistance[I]) < GREAT)
        {
            signedDistance[I] = frontSurface.signedDistance(C[I]);
        }
    }
    
    // cell-corner distances
    auto& pointSignedDistance = lookupPointSignedDistance();
    const auto& P = mesh().points();

    forAll(pointSignedDistance, I)
    {
        if (mag(pointSignedDistance[I]) < GREAT)
        {
            pointSignedDistance[I] = frontSurface.signedDistance(P[I]);
        }
    }
}

void lentSubalgorithmTest::computeFrontSignedDistances()
{
    auto& signedDistance = lookupSignedDistance();
    auto& pointSignedDistance = lookupPointSignedDistance();

    const auto& searchDistanceSqr = lookupSearchDistanceSqr();
    const auto& pointSearchDistanceSqr = lookupPointSearchDistanceSqr();

    lent_.calcSignedDistances
            (
                signedDistance,
                pointSignedDistance,
                searchDistanceSqr,
                pointSearchDistanceSqr,
                front_
            ); 
}

void lentSubalgorithmTest::setupFrontFromSurface(const bool correct)
{
    surfaceTmp_.ref().randomize();
    surfaceTmp_.ref().writeParameters(assembleFilePath() + "front/surfaceParameters.dat");

    computeExactSignedDistances();

    lent_.reconstructFront(front_, lookupSignedDistance(), lookupPointSignedDistance());

    if (correct)
    {
        surfaceTmp_.ref().moveFrontToSurface(front_);
        surfaceTmp_.ref().makeNormalOrientationConsistent(front_);
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentSubalgorithmTest::lentSubalgorithmTest(const fvMesh& mesh, triSurfaceFront& front)
:
    mesh_{mesh},
    front_{front},
    lent_{front, mesh},
    surfaceTmp_{},
    separator_{","},
    algorithmRuntime_{"t_exec"}
{
    surfaceTmp_ = tmp<analyticalSurface>
                  {
                     analyticalSurface::New(lentDict().subDict("frontSurface"))
                  };

    nRandomRuns_ = testDict().get<label>("nRandomRuns");
    nPerturbedRuns_ = testDict().get<label>("nPerturbedRuns");
    writeFields_ = testDict().get<Switch>("writeFields");

    // Values below 1 for the number of random runs and for the number
    // of perturbed runs make no sense.
    if (nRandomRuns_ < 1 || nPerturbedRuns_ < 1)
    {
        FatalErrorIn
        (
            "lentSubalgorithmTest::lentSubalgorithmTest(const fvMesh& mesh, triSurfaceFront& front)"
        )   << "Invalid entry for nRandomRuns or nPerturbedRuns: "
            << "Use n >= 1" 
            << abort(FatalError);
    }

    // Initialize search distances
    auto& searchDistanceSqr = lookupSearchDistanceSqr();
    auto& pointSearchDistanceSqr = lookupPointSearchDistanceSqr();

    lent_.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void lentSubalgorithmTest::runAllTests(const word& fileName)
{
    Info << "Running Tests..." << endl;

    for (label I = 0; I < nRandomRuns_; ++I)
    {
        Info << "\n\n------------------------------------------------------"
             << "\n ---> Running random setup...\n";
        randomSetup();

        for (label K = 0; K < nPerturbedRuns_; ++K)
        {
            Info << "\n---> Perturbing input fields for model...\n";
            
            // Use the Time class to write the fields and front of
            // each run if desired
            setNextRun();
        
            perturbInputFields();

            // add time measurement
            auto start = clock::now();

            Info << "\n---> Running model...\n";
            computeApproximatedFields();

            auto end = clock::now();
            scalar deltaT = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            addMeasure(algorithmRuntime_, deltaT);
            
            Info << "\n---> Evaluating test metrics...\n";
            evaluateMetrics();

            writeFields();
            writeMetrics(fileName, I*nPerturbedRuns_ + K);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
