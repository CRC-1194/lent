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

#include "analyticalSurface.H"
#include "lentMethod.H"
#include "triSurfaceFront.H"

#include <string>
#include <unordered_map>
#include <vector>

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
    template<class Type>
    using metricDataMap = std::unordered_map<std::string, std::vector<Type>>;

    // Private data
    const fvMesh& mesh_;
    triSurfaceFront& front_;
    lentMethod lent_;
    tmp<analyticalSurface> surfaceTmp_;

    label nRandomRuns_;
    label nPerturbedRuns_;
    Switch writeFields_;
    std::string separator_;
    std::string algorithmRuntime_;

    metricDataMap<scalar> scalarMetrics_;
    metricDataMap<vector> vectorMetrics_;

    // Private Member Functions
    std::string assembleFilePath() const;
    template<typename Type>
    std::string metricHeader(const metricDataMap<Type>&) const;
    std::string printFoamVector(const vector& v) const;
    void writeMetrics(const word& fileName, const unsigned long runNumber) const;

    void setNextRun() const;
    void writeFields() const;

    volScalarField& lookupSearchDistanceSqr() const
    {
        return lookupField<volScalarField>("searchDistanceSqr");
    }

    pointScalarField& lookupPointSearchDistanceSqr() const
    {
        return lookupField<pointScalarField>("pointSearchDistanceSqr");
    }


protected:
    void addMeasure(const std::string& metricName, const scalar& value);
    void addMeasure(const std::string& metricName, const vector& value);

    const fvMesh& mesh() const
    {
        return mesh_;
    }

    triSurfaceFront& frontRef()
    {
        return front_;
    }

    const triSurfaceFront& frontRef() const
    {
        return front_;
    }

    auto& lent()
    {
        return lent_;
    }

    const dictionary& lentDict() const
    {
        return lent_.dict();
    }

    const dictionary& testDict() const
    {
        return lentDict().subDict("testConfiguration");
    } 

    label nRandomRuns() const
    {
        return nRandomRuns_;
    }

    label nPerturbedRuns() const
    {
        return nPerturbedRuns_;
    }

    const analyticalSurface& surfaceRef() const
    {
        return surfaceTmp_.ref();
    }

    template<class fieldType>
    fieldType& lookupField(const word fieldName) const
    {
        return const_cast<fieldType&>(mesh_.lookupObject<fieldType>(fieldName));
    }

    volScalarField& lookupSignedDistance() const
    {
        return lookupField<volScalarField>("signedDistance");
    }
    
    pointScalarField& lookupPointSignedDistance() const
    {
        return lookupField<pointScalarField>("pointSignedDistance");
    }

    void computeExactSignedDistances() const;
    void computeExactSignedDistancesNarrowBand() const;
    void computeFrontSignedDistances();
    void setupFrontFromSurface(const bool correct);


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
    void runAllTests(const word& fileName);
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
