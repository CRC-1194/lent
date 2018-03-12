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
#include <sstream>
#include <string>

namespace Foam {
namespace FrontTracking {

    using clock = std::chrono::high_resolution_clock;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
lentSubalgorithmTest::fileFormat lentSubalgorithmTest::detectFileFormat(const word& fileName) const
{
    if (fileName.rfind(".hdf5") != std::string::npos)
    {
        return fileFormat::hdf5;
    }
    else if (fileName.rfind(".csv") != std::string::npos)
    {
        return fileFormat::csv;
    }
    else
    {
        return fileFormat::fallBack;
    }
}

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
    auto& runTime = const_cast<Time&>(mesh_.time());
    runTime.write();
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
    const auto& C = mesh_.C();
    auto& signedDistance = lookupSignedDistance();

    forAll(signedDistance, I)
    {
        signedDistance[I] = frontSurface.signedDistance(C[I]);
    }
    
    // cell-corner distances
    auto& pointSignedDistance = lookupPointSignedDistance();
    const auto& pMesh = pointSignedDistance.mesh();
    const auto& corners = pMesh().points();

    forAll(pointSignedDistance, I)
    {
        pointSignedDistance[I] = frontSurface.signedDistance(corners[I]);
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
    // re-read analyticalSurface in order to randomize it
    surfaceTmp_ = tmp<analyticalSurface>
                  {
                        analyticalSurface::New(lentDict().subDict("frontSurface"))
                  };

    computeExactSignedDistances();

    lent_.reconstructFront(front_, lookupSignedDistance(), lookupPointSignedDistance());

    if (correct)
    {
        surfaceTmp_.ref().moveFrontToSurface(front_);   
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

    nRandomRuns_ = lentDict().lookupOrDefault<label>("nRandomRuns", 1);
    nPerturbedRuns_ = lentDict().lookupOrDefault<label>("nPerturbedRuns", 1);

    // Values below 1 for the number of random runs and for the number
    // of perturbed runs make no sense.
    // Rather to throw an error, set the values to 1 if they are 0 or smaller
    if (nRandomRuns_ < 1)
    {
       nRandomRuns_ = 1;
    }

    if (nPerturbedRuns_ < 1)
    {
       nPerturbedRuns_ = 1;
    } 

    // Initialize search distances
    auto& searchDistanceSqr = lookupSearchDistanceSqr();
    auto& pointSearchDistanceSqr = lookupPointSearchDistanceSqr();

    lent_.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);
}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void lentSubalgorithmTest::runAllTests()
{
    Info << "Running Tests..." << endl;

    for (label I = 0; I < nRandomRuns_; ++I)
    {
        // Use the Time class to write the fields and front of
        // each random run if desired
        setNextRun();
        
        randomSetup();

        for (label K = 0; K < nPerturbedRuns_; ++K)
        {
            perturbInputFields();

            // add time measurement
            auto start = clock::now();

            computeApproximatedFields();

            auto end = clock::now();
            scalar deltaT = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            addMeasure(algorithmRuntime_, deltaT);
            
            evaluateMetrics();
        }

        writeFields();
    }
}

void lentSubalgorithmTest::writeResults(const word& fileName) const
{
    Info << "Writing results..." << endl;

    auto fileType = detectFileFormat(fileName);
    
    if (fileType == fileFormat::hdf5)
    {
        writeResultsHDF5(fileName);
    }
    else if (fileType == fileFormat::csv)
    {
        writeResultsCSV(fileName);
    }
    else if (fileType == fileFormat::fallBack)
    {
        // Use csv as default file format in case there is no specialized
        // writer. The computations have been performed anyway, so do not
        // throw them away (TT)
        word fallBackFileName{fileName + ".csv"};

        Info << "NOTE: no writer found for " << fileName << ". Using "
             << fallBackFileName << " instead.\n";

        writeResultsCSV(fallBackFileName);
    }
}

void lentSubalgorithmTest::writeResultsHDF5(const word& fileName) const
{
    notImplemented("writeResultsHDF5(...)");
}

void lentSubalgorithmTest::writeResultsCSV(const word& fileName) const
{
    std::fstream dataFile(assembleFilePath() + fileName, std::ios_base::out);

    dataFile << metricHeader(scalarMetrics_)
             << metricHeader(vectorMetrics_) << std::endl;

    // TODO: assumption: for every run every metric is evaluated,
    //  so that the results will be a full table without empty fields
    label nRuns = scalarMetrics_.at(algorithmRuntime_).size();

    for (int index = 0; index < nRuns; ++index)
    {
        for (auto const& metricField: scalarMetrics_)
        {
            dataFile << metricField.second[index] << separator_;
        } 

        for (auto const& metricField: vectorMetrics_)
        {
            dataFile << printFoamVector(metricField.second[index]) << separator_;
        } 
        
        dataFile << std::endl;
    }

    dataFile.close();
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
