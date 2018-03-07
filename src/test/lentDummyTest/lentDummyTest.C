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

#include "errorMetrics.H"

#include <chrono>
#include <thread>

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void lentDummyTest::randomSetup()
{
}

void lentDummyTest::computeApproximatedFields()
{
    auto& computedField = computedFieldTmp_.ref();

    forAll(computedField, I)
    {
        computedField[I] = I;
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(100));
}

void lentDummyTest::evaluateMetrics()
{
    auto deltaField = computedFieldTmp_.ref() - referenceFieldTmp_.ref();

    errorMetrics eval{deltaField.ref()};

    addMeasure("L1_norm", eval.arithmeticMeanError());
    addMeasure("L2_norm", eval.quadraticMeanError());
    addMeasure("L_inf_norm", eval.maximumError());

    // Write fields
    referenceFieldTmp_.ref().write();
    computedFieldTmp_.ref().write();
    deltaField.ref().write();
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
            "refField",
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
            "computedField",
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
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
