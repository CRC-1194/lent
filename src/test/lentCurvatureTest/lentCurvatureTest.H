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
    Foam::lentCurvatureTest

Description

SourceFiles
    lentCurvatureTestI.H
    lentCurvatureTest.C
    lentCurvatureTestIO.C

\*---------------------------------------------------------------------------*/

#ifndef lentCurvatureTest_H
#define lentCurvatureTest_H

#include "frontExactSurfaceCurvatureModel.H"
#include "lentSubalgorithmTest.H"
#include "noiseGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class lentCurvatureTest Declaration
\*---------------------------------------------------------------------------*/

class lentCurvatureTest
:
    public lentSubalgorithmTest
{
    // Private data
    tmp<volScalarField> filterFieldTmp_;
    tmp<surfaceScalarField> faceFilterFieldTmp_;
    tmp<volScalarField> relativeDeltaFieldTmp_;
    tmp<surfaceScalarField> relativeFaceDeltaFieldTmp_;

    noiseGenerator noiseGen_;
    tmp<frontExactSurfaceCurvatureModel> exactCurvatureModelTmp_;
    tmp<frontCurvatureModel> numericalCurvatureModelTmp_;

    Switch correctFront_;
    Switch useFrontSignedDistance_;
    vector frontNoise_;
    scalar distanceNoise_;


    // Private Member Functions
    virtual void randomSetup();
    virtual void perturbInputFields();
    virtual void computeApproximatedFields();
    virtual void evaluateMetrics();


public:

    // Constructors
    lentCurvatureTest(const fvMesh& mesh, triSurfaceFront& front);


    //- Destructor
    virtual ~lentCurvatureTest() = default;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
