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
    Foam::analyticalSurfaceNormalConsistency

SourceFiles
    analyticalSurfaceNormalConsistency.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Uses an analytical surface description to ensure consistent normal
    orientation.
    Intended for testing purposes only.

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/

#include "analyticalSurfaceNormalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalSurfaceNormalConsistency, 0);
    addToRunTimeSelectionTable(normalConsistency, analyticalSurfaceNormalConsistency, Dictionary);

// * * * * * * * * * * * * * Private  member functions * * * * * * * * * * * //
void analyticalSurfaceNormalConsistency::runNormalConsistencyAlgorithm(
    triSurfaceFront& front,
    const volScalarField&,
    const pointScalarField&
) const
{
    makeFrontNormalsConsistent(front);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalSurfaceNormalConsistency::analyticalSurfaceNormalConsistency(const dictionary& configDict)
    :
        normalConsistency{configDict},
        reverseNormal_{configDict.lookupOrDefault<bool>("reverseNormal", false)},
        surfaceTmp_{analyticalSurface::New(configDict.subDict("frontSurface"))}
{
}
    
analyticalSurfaceNormalConsistency::analyticalSurfaceNormalConsistency(const analyticalSurface& surface, const bool reverseNormal)
:
    normalConsistency{dictionary{}},
    reverseNormal_{reverseNormal},
    surfaceTmp_{surface}
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool analyticalSurfaceNormalConsistency::makeFrontNormalsConsistent(triSurfaceFront& front) const
{
    bool normalsWereConsistent = true;

    List<labelledTri>& triangles = static_cast<List<labelledTri>& > (front);
    const auto& triangleNormals = front.faceNormals();
    const auto& faceCentres = front.Cf();
    const auto& surface = surfaceTmp_.operator()();

    scalar normalSign = 1.0;

    if (reverseNormal_)
    {
        normalSign = -1.0;
    }

    // For all faces
    forAll (triangles, triangleI)
    {
        auto normal = normalSign*surface.normalToPoint(faceCentres[triangleI]);
        if ((normal & triangleNormals[triangleI]) < 0.0)
        {
            triangles[triangleI].flip();
            normalsWereConsistent = false;
        }
    }

    return normalsWereConsistent;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
