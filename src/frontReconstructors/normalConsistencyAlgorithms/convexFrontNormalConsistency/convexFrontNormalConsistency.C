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
    Foam::convexFrontNormalConsistency

SourceFiles
    convexFrontNormalConsistency.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Ensures normal consistency by computing a geometric centre.
    ONLY applicable to CONVEX fronts.

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

#include "convexFrontNormalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(convexFrontNormalConsistency, 0);
    addToRunTimeSelectionTable(normalConsistency, convexFrontNormalConsistency, Dictionary);

// * * * * * * * * * * * * * Private  member functions * * * * * * * * * * * //
void convexFrontNormalConsistency::runNormalConsistencyAlgorithm(
    triSurfaceFront& front,
    const volScalarField&,
    const pointScalarField&
) const
{
    List<labelledTri>& triangles = static_cast<List<labelledTri>& > (front);
    const auto& triangleNormals = front.faceNormals();
    const auto& points = front.points();
    const auto& faceCentres = front.Cf();

    vector geometricCentre{0.0, 0.0, 0.0};

    forAll(points, I)
    {
        geometricCentre += points[I];
    }

    geometricCentre /= points.size();

    // For all faces
    forAll (triangles, triangleI)
    {
        auto centreToPoint = faceCentres[triangleI] - geometricCentre;

        if (sign(centreToPoint & triangleNormals[triangleI]) != orientationSign_)
        {
            triangles[triangleI].flip();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

convexFrontNormalConsistency::convexFrontNormalConsistency(const dictionary& configDict)
    :
        normalConsistency(configDict),
        orientation_{configDict.get<word>("normalOrientation")},
        orientationSign_{}
{
    if (orientation_ == "outside")
    {
        orientationSign_ = 1.0;
    }
    else if (orientation_ == "inside")
    {
        orientationSign_ = -1.0;
    }
    else
    {
        FatalErrorIn (
            "convexFrontNormalConsistency::convexFrontNormalConsistency(const dictionary& configDict)"
        )   << "Unknown normalOrientation "
            << orientation_ << nl << nl
            << "Valid orientations are : " << endl
            << "outside or inside"
            << exit(FatalError);
    }
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
