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
    Foam::diffuseInterfaceProperties

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

Description
    
    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

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
