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
