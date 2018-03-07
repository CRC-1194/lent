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

Class
    Foam::lentSubalgorithmTest

Description

SourceFiles
    lentSubalgorithmTestI.H
    lentSubalgorithmTest.C
    lentSubalgorithmTestIO.C

\*---------------------------------------------------------------------------*/

#ifndef lentSubalgorithmTest_H
#define lentSubalgorithmTest_H

#include "fvMesh.H"
#include "messageStream.H"
#include "vector.H"
#include "word.H"

#include "lentMethod.H"
#include "triSurfaceFront.H"

#include <string>
#include <unordered_map>
#include <vector>

// Forward declarations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class lentSubalgorithmTest Declaration
\*---------------------------------------------------------------------------*/

class lentSubalgorithmTest
{

private:
    // Aliases
    enum class fileFormat {csv, hdf5, fallBack};

    template<class Type>
    using metricDataMap = std::unordered_map<std::string, std::vector<Type>>;

    // Private data
    const fvMesh& mesh_;
    triSurfaceFront& front_;
    lentMethod lent_;

    std::string separator_;
    std::string algorithmRuntime_;

    metricDataMap<scalar> scalarMetrics_;
    metricDataMap<vector> vectorMetrics_;

    // Private Member Functions
    fileFormat detectFileFormat(const word& fileName) const;
    std::string assembleFilePath() const;
    template<typename Type>
    std::string metricHeader(const metricDataMap<Type>&) const;
    std::string printFoamVector(const vector& v) const;


protected:
    void addMeasure(const std::string& metricName, const scalar& value);
    void addMeasure(const std::string& metricName, const vector& value);
    auto& lent()
    {
        return lent_;
    }

    virtual void randomSetup() = 0;
    virtual void perturbInputFields() = 0;
    virtual void computeApproximatedFields() = 0;
    virtual void evaluateMetrics() = 0;


public:

    // Constructors
    lentSubalgorithmTest(const fvMesh& mesh, triSurfaceFront& front);


    //- Destructor
    virtual ~lentSubalgorithmTest() = default;


    // Member Functions
    void runAllTests();
    void writeResults(const word& fileName) const;
    void writeResultsHDF5(const word& fileName) const;
    void writeResultsCSV(const word& fileName) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lentSubalgorithmTestI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
