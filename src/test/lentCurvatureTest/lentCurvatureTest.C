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

#include "lentCurvatureTest.H"

#include "fvCFD.H"
#include "volMesh.H"

#include "errorMetrics.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void lentCurvatureTest::randomSetup()
{
    setupFrontFromSurface(correctFront_);

    // Remember: computing the signed distance also updates front-mesh
    // connectivity
    computeFrontSignedDistances();

    if (!useFrontSignedDistance_)
    {
        computeExactSignedDistancesNarrowBand();
    }

    // The markerfield is required to setup the filter field
    auto& markerField = lookupField<volScalarField>("alpha.water");
    lent().calcMarkerField(markerField);

    auto& filterField = filterFieldTmp_.ref();
    dimensionedScalar dSMALL("SMALL", pow(dimLength,-1), SMALL);
    filterField = pos(mag(fvc::grad(markerField)) - 1e2*dSMALL);

    auto& faceFilterField = faceFilterFieldTmp_.ref();
    faceFilterField = mag(fvc::snGrad(markerField))
                        * dimensionedScalar{"one", dimLength, 1.0};

    forAll(faceFilterField, I)
    {
        if (faceFilterField[I] > SMALL)
        {
            faceFilterField[I] = 1.0;
        }
    }

    // Boundary
    // TODO: Think about reasonable filtering approach once we start to treat
    // contact line problems.
    // For now, disable boundary faces since their snGrad depends on the
    // boundary condition of the markerfield
    forAll(faceFilterField.boundaryFieldRef(), I)
    {
        auto& bField = faceFilterField.boundaryFieldRef()[I];

        forAll(bField, K)
        {
            if (bField[K] > SMALL)
            {
                bField[K] = 0.0;
            }
        }
    }
}

void lentCurvatureTest::perturbInputFields()
{
    // For approaches working directly on the front, disturb front vertices.
    // Perturbation only makes sense when applied to the exact vertex
    // positions
    if (correctFront_ && mag(frontNoise_) > SMALL)
    {
        const auto& surface = surfaceRef();
        auto& front = frontRef();

        surface.moveFrontToSurface(front);
        auto& vertices = const_cast<pointField&>(front.points());
        noiseGen_.addNoiseTo<vector,List>(vertices, frontNoise_);
        
       // Clear buffered point dependent data
       front.clearGeom();
    }

    // Analogue to the front based approaches, perform the same for the signed
    // distance field
    if (!useFrontSignedDistance_ && distanceNoise_ > SMALL)
    {
        auto& signedDistance = lookupSignedDistance();

        computeExactSignedDistances();

        noiseGen_.addNoiseTo<scalar,List>(signedDistance, distanceNoise_);
    }
}

void lentCurvatureTest::computeApproximatedFields()
{
    const auto& numericalCurvatureModel = numericalCurvatureModelTmp_.ref();
    numericalCurvatureModel.cellCurvature(mesh(), frontRef());
}

void lentCurvatureTest::evaluateMetrics()
{
    // Get fields and compute error field
    auto& exactCurvatureModel = exactCurvatureModelTmp_.ref();
    auto& numericalCurvatureModel = numericalCurvatureModelTmp_.ref();
    auto exactCurvatureFieldPtr = exactCurvatureModel.cellCurvature(mesh(), frontRef());
    auto& exactCurvatureField = *exactCurvatureFieldPtr;
    auto numericalCurvatureFieldPtr = numericalCurvatureModel.cellCurvature(mesh(), frontRef());
    auto& numericalCurvatureField = *numericalCurvatureFieldPtr;

    const auto& filterField = filterFieldTmp_.ref();
    exactCurvatureField *= filterField;
    numericalCurvatureField *= filterField;

    auto& relativeDeltaField = relativeDeltaFieldTmp_.ref();
    dimensionedScalar dSMALL("SMALL", pow(dimLength,-1), SMALL);

    relativeDeltaField = mag(numericalCurvatureField - exactCurvatureField)/(mag(exactCurvatureField) + dSMALL);

    // Rename curvature fields so they can be distinguished
    exactCurvatureField.rename("exact_cell_curvature");
    numericalCurvatureField.rename("numerical_cell_curvature");

    // Acutually compute metrics
    auto minNumerical = min(numericalCurvatureField);
    auto maxNumerical = max(numericalCurvatureField);
    auto minExact = min(exactCurvatureField);
    auto maxExact = max(exactCurvatureField);

    errorMetrics eval{relativeDeltaField};

    addMeasure("L1_norm", eval.arithmeticMeanError());
    addMeasure("L2_norm", eval.quadraticMeanError());
    addMeasure("L_inf_norm", eval.maximumError());
    addMeasure("min_numerical", minNumerical.value());
    addMeasure("max_numerical", maxNumerical.value());
    addMeasure("min_exact", minExact.value());
    addMeasure("max_exact", maxExact.value());

    // Evaluation for curvature at faces -> location where the curvature
    // is needed for discretization of the surface tension force
    auto exactFaceCurvaturePtr = exactCurvatureModel.faceCurvature(mesh(), frontRef());
    const auto& exactFaceCurvature = *exactFaceCurvaturePtr;
    auto numericalFaceCurvaturePtr = numericalCurvatureModel.faceCurvature(mesh(), frontRef());
    const auto& numericalFaceCurvature = *numericalFaceCurvaturePtr;
    auto& relativeFaceDeltaField = relativeFaceDeltaFieldTmp_.ref();
    
    relativeFaceDeltaField = mag(numericalFaceCurvature - exactFaceCurvature)/(mag(exactFaceCurvature) + dSMALL);
    relativeFaceDeltaField *= faceFilterFieldTmp_.ref();
    
    errorMetrics faceEval{relativeFaceDeltaField};
    addMeasure("face L1_norm", faceEval.arithmeticMeanError());
    addMeasure("face L2_norm", faceEval.quadraticMeanError());
    addMeasure("face L_inf_norm", faceEval.maximumError());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentCurvatureTest::lentCurvatureTest(const fvMesh& mesh, triSurfaceFront& front)
:
    lentSubalgorithmTest{mesh, front},
    filterFieldTmp_{},
    relativeDeltaFieldTmp_{},
    noiseGen_{},
    exactCurvatureModelTmp_{},
    numericalCurvatureModelTmp_{}
{
    correctFront_ = Switch{testDict().lookup("correctFront")};
    useFrontSignedDistance_ = Switch{testDict().lookup("useFrontSignedDistance")};
    frontNoise_ = testDict().lookup("frontNoise");
    distanceNoise_ = readScalar(testDict().lookup("distanceNoise"));

    exactCurvatureModelTmp_ = tmp<analyticalSurfaceCurvatureModel>{
        new analyticalSurfaceCurvatureModel(lentDict().subDict("exactCurvatureModel"), surfaceRef())
    };

    numericalCurvatureModelTmp_ = tmp<curvatureModel>{
        curvatureModel::New(lentDict().subDict("numericalCurvatureModel"))
    };

    // Initialize Fields
    filterFieldTmp_ = tmp<volScalarField>{new volScalarField 
    (
        IOobject
        (
            "filter_field",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(
            "zero",
            dimless,
            0.0
        )
    )
    };

    faceFilterFieldTmp_ = tmp<surfaceScalarField>{new surfaceScalarField
    (
        IOobject(
            "face_filter_field", 
            mesh.time().timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ), 
        mesh, 
        dimensionedScalar(
            "zero", 
            dimless, 
            0.0
        )
    )
    };

    relativeDeltaFieldTmp_ = tmp<volScalarField>{new volScalarField 
    (
        IOobject
        (
            "relative_curvature_error",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(
            "zero",
            dimless,
            0.0
        )
    )
    };

    relativeFaceDeltaFieldTmp_ = tmp<surfaceScalarField>{new surfaceScalarField
    (
        IOobject(
            "face_relative_curvature_error", 
            mesh.time().timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ), 
        mesh, 
        dimensionedScalar(
            "zero", 
            dimless, 
            0.0
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
