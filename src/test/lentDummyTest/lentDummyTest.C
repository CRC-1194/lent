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

#include "lentDummyTest.H"

#include "analyticalSurface.H"
#include "errorMetrics.H"

#include <chrono>
#include <thread>

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void lentDummyTest::randomSetup()
{
    Info << "\n\n----------------------------------------------\n"
         << "Setup...\n";

    setupFrontFromSurface(true);
    computeFrontSignedDistances();
}

void lentDummyTest::perturbInputFields()
{
    Info << "Perturbing fields...\n";

    auto& referenceField = referenceFieldTmp_.ref();
    auto& refVectorField = refVectorFieldTmp_.ref();

    referenceField = dimensionedScalar{"zero", dimless, 0.0};
    refVectorField = dimensionedVector{"one", dimless, vector{1,1,0}};

    noiseGen_.addNoiseTo<scalar,List>(referenceField, 1.0);
    noiseGen_.addNoiseTo<vector,List>(refVectorField, vector{1,0,0});
}

void lentDummyTest::computeApproximatedFields()
{
    Info << "Approximating fields...\n";

    auto& computedField = computedFieldTmp_.ref();
    auto& compVectorField = compVectorFieldTmp_.ref();

    forAll(computedField, I)
    {
        computedField[I] = 10.0;
        compVectorField[I] = vector{10, 1, 0};
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(50));
}

void lentDummyTest::evaluateMetrics()
{
    Info << "Evaluating metrics...\n";

    auto deltaField = computedFieldTmp_.ref() - referenceFieldTmp_.ref();
    auto deltaVetcorFieldTmp = compVectorFieldTmp_.ref() - refVectorFieldTmp_.ref();

    errorMetrics eval{deltaField.ref()};

    addMeasure("L1_norm", eval.arithmeticMeanError());
    addMeasure("L2_norm", eval.quadraticMeanError());
    addMeasure("L_inf_norm", eval.maximumError());

    // Write delta fields
    deltaField.ref().write();
    deltaVetcorFieldTmp.ref().write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentDummyTest::lentDummyTest(const fvMesh& mesh, triSurfaceFront& front)
:
    lentSubalgorithmTest{mesh, front}
{
    referenceFieldTmp_ = tmp<volScalarField>{new volScalarField 
    (
        IOobject
        (
            "refScalarField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(
            "zero",
            dimless,
            0
        )
    )
    };

    computedFieldTmp_ = tmp<volScalarField>{new volScalarField 
    (
        IOobject
        (
            "computedScalarField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(
            "zero",
            dimless,
            1
        )
    )
    };

    refVectorFieldTmp_ = tmp<volVectorField>{new volVectorField 
    (
        IOobject
        (
            "refVectorField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(
            "zero",
            dimless,
            vector{1,1,0}
        )
    )
    };

    compVectorFieldTmp_ = tmp<volVectorField>{new volVectorField 
    (
        IOobject
        (
            "compVectorField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(
            "zero",
            dimless,
            vector{1,1,0}
        )
    )
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
