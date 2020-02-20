/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::lentGlobalMarkerFieldTest

SourceFiles
    lentGlobalMarkerFieldTest.C

Author
    Tobias Tolle   tolle@mma.tu-darmstadt.de

Description
    This class defines markerfield tests on a full domain scale.
    The following properties are tested:
        * boundedness of the marker field: alpha in [0;1]
        * relative global volume error based on the front volume

    For now, this test does not make use of exact signed distances and assumes
    that the alpha value of the dispersed phase is 0.0.


    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#include "lentGlobalMarkerFieldTest.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
bool lentGlobalMarkerFieldTest::markerFieldIsBounded(const volScalarField& markerField) const
{
    for (const auto& alpha : markerField)
    {
        if (alpha < 0.0 || alpha > 1.0)
        {
            return false;
        }
    }

    return true;
}

scalar lentGlobalMarkerFieldTest::computeMarkerFieldVolume(const volScalarField& markerField) const
{
    scalar volume{0.0};

    const auto& cellVolumes = mesh().V();

    forAll(markerField, I)
    {
        volume += cellVolumes[I] * (1.0 - markerField[I]);
    }

    return volume;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void lentGlobalMarkerFieldTest::randomSetup()
{
    setupFrontFromSurface(correctFront_);

    computeFrontSignedDistances();
}

void lentGlobalMarkerFieldTest::perturbInputFields()
{
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
       
        computeFrontSignedDistances();
    }
}

void lentGlobalMarkerFieldTest::computeApproximatedFields()
{
    auto& markerField = lookupField<volScalarField>(phaseIndicatorFieldName_);
    lent().calcMarkerField(markerField);
}

void lentGlobalMarkerFieldTest::evaluateMetrics()
{
    const auto frontVolume = frontRef().convexFrontVolume();

    const auto& markerField = lookupField<volScalarField>(phaseIndicatorFieldName_);
    const auto markerFieldVolume = computeMarkerFieldVolume(markerField);
    const auto fieldIsBounded = markerFieldIsBounded(markerField);

    const auto relativeVolumeError = mag(markerFieldVolume - frontVolume) / (frontVolume + SMALL);

    addMeasure("is_bounded", scalar(fieldIsBounded));
    addMeasure("relative_volume_error", relativeVolumeError);
    addMeasure("front_volume", frontVolume);
    addMeasure("marker_field_volume", markerFieldVolume);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentGlobalMarkerFieldTest::lentGlobalMarkerFieldTest(const fvMesh& mesh, triSurfaceFront& front)
:
    lentSubalgorithmTest{mesh, front},
    correctFront_{},
    frontNoise_{},
    phaseIndicatorFieldName_{}
{
    correctFront_ = testDict().get<Switch>("correctFront");
    frontNoise_ = testDict().get<vector>("frontNoise");
    phaseIndicatorFieldName_ = testDict().get<word>("phaseIndicatorFieldName");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
