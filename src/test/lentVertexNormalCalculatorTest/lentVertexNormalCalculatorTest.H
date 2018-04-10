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
    Foam::lentVertexNormalCalculatorTest

Description

SourceFiles
    lentVertexNormalCalculatorTestI.H
    lentVertexNormalCalculatorTest.C
    lentVertexNormalCalculatorTestIO.C

\*---------------------------------------------------------------------------*/

#ifndef lentVertexNormalCalculatorTest_H
#define lentVertexNormalCalculatorTest_H

#include "frontVertexNormalCalculator.H"
#include "lentSubalgorithmTest.H"
#include "noiseGenerator.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class lentVertexNormalCalculatorTest Declaration
\*---------------------------------------------------------------------------*/

class lentVertexNormalCalculatorTest
:
    public lentSubalgorithmTest
{

private:

    // Private data
    noiseGenerator noiseGen_{};

    tmp<frontVertexNormalCalculator> normalCalculatorTmp_;

    Switch correctFront_;
    scalar magFrontNoise_;

    // Data for metric evaluation
    tmp<triSurfaceFrontPointVectorField> approximateNormalsTmp_;
    std::vector<scalar> normalDeviationField_;


    // Private Member Functions
    // TODO: there are overlaps between this test class and the front smoothing
    // test class. Remove duplication when refactoring
    tensor twoDProjector() const;
    void addFrontNormalNoise();


protected:
    virtual void randomSetup();
    virtual void perturbInputFields();
    virtual void computeApproximatedFields();
    virtual void evaluateMetrics();

public:

    // Constructors
    lentVertexNormalCalculatorTest(const fvMesh& mesh, triSurfaceFront& front);


    //- Destructor
    ~lentVertexNormalCalculatorTest() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
