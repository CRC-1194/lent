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
    Foam::lentFrontSmoothingTest

Description

SourceFiles
    lentFrontSmoothingTestI.H
    lentFrontSmoothingTest.C
    lentFrontSmoothingTestIO.C

\*---------------------------------------------------------------------------*/

#ifndef lentFrontSmoothingTest_H
#define lentFrontSmoothingTest_H

#include <memory>

#include "frontSmoother.H"
#include "lentSubalgorithmTest.H"
#include "noiseGenerator.H"
#include "analyticalSurfaceNormalConsistency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class lentFrontSmoothingTest Declaration
\*---------------------------------------------------------------------------*/

class lentFrontSmoothingTest
:
    public lentSubalgorithmTest
{
    // Private data
    noiseGenerator noiseGen_{};

    std::unique_ptr<frontSmoother> frontSmootherPtr_;
    std::unique_ptr<analyticalSurfaceNormalConsistency> normalConsistencyPtr_;

    Switch correctFront_;
    scalar magFrontNoise_;
    
    // Data for metric evaluation
    scalar originalVolume_;
    point originalCentreOfGravity_;
    tmp<triSurfaceFrontPointScalarField> surfaceDeviationTmp_;


    // Private Member Functions
    bool normalsAreConsistent();
    bool boundaryPointsRemainOnBoundary() const;

    // TODO: volume computation should be moved to triSurfaceFront,
    // named like "convexFrontVolume" so that it is clear that the
    // approach only works for convex surfaces
    bool is2DCase() const;
    point edgeLoopCentre(const labelList& loop, const triSurfaceFront& front) const;
    scalar tetVolume(const point& base, const point& v0, const point& v1, const point& v2) const;
    tensor twoDProjector() const;
    void addFrontNormalNoise();

    scalar frontVolume() const;
    point computeFrontCentreOfGravity() const;
    void computeSurfaceDeviation();
    
    virtual void randomSetup();
    virtual void perturbInputFields();
    virtual void computeApproximatedFields();
    virtual void evaluateMetrics();


public:

    // Constructors
    lentFrontSmoothingTest(const fvMesh& mesh, triSurfaceFront& front);


    //- Destructor
    virtual ~lentFrontSmoothingTest() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
