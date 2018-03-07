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

#include <iostream>

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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void lentSubalgorithmTest::addMeasure(const std::string& metricName, const scalar& value)
{
    scalarMetrics_[metricName].push_back(value);
}

void lentSubalgorithmTest::addMeasure(const std::string& metricName, const vector& value)
{
    vectorMetrics_[metricName].push_back(value);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentSubalgorithmTest::lentSubalgorithmTest(const fvMesh& mesh, triSurfaceFront& front)
:
    mesh_{mesh},
    front_{front},
    lent_{front, mesh},
    separator_{","},
    algorithmRuntime_{"t_exec"}
{}



// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void lentSubalgorithmTest::runAllTests()
{
    Info << "Running Tests..." << endl;

    label nRandomRuns = 2;
    label nPerturbedRuns = 2;

    for (label I = 0; I < nRandomRuns; ++I)
    {
        randomSetup();

        for (label K = 0; K < nPerturbedRuns; ++K)
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
    std::string filePath = assembleFilePath() + fileName;

    std::cout << "File path = " << filePath << std::endl;

    std::fstream dataFile(filePath, std::ios_base::out);

    dataFile << metricHeader(scalarMetrics_)
             << metricHeader(vectorMetrics_) << std::endl;

    // TODO: assumption: for every run every mertic is evaluated,
    //  so that the results will be a full table without empty fields
    label nRuns = scalarMetrics_.at(algorithmRuntime_).size();

    // Skip first entry since it is not obtained by evaluation but by
    // default initialization
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
