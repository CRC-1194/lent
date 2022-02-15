/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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

\*---------------------------------------------------------------------------*/

#ifndef lentGlobalMarkerFieldTest_H
#define lentGlobalMarkerFieldTest_H

#include "lentSubalgorithmTest.H"
#include "noiseGenerator.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class lentGlobalMarkerFieldTest Declaration
\*---------------------------------------------------------------------------*/

class lentGlobalMarkerFieldTest
:
    public lentSubalgorithmTest
{

private:

    // Private data
    noiseGenerator noiseGen_;
    
    Switch correctFront_;
    vector frontNoise_;
    word phaseIndicatorFieldName_;

    // Private Member Functions
    bool markerFieldIsBounded(const volScalarField& markerField) const;
    scalar computeMarkerFieldVolume(const volScalarField& markerField) const;
    
protected:

    virtual void randomSetup();
    virtual void perturbInputFields();
    virtual void computeApproximatedFields();
    virtual void evaluateMetrics();


public:

    // Constructors
    lentGlobalMarkerFieldTest(const fvMesh& mesh, triSurfaceFront& front);


    //- Destructor
    virtual ~lentGlobalMarkerFieldTest() = default;


    // Member Functions
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
