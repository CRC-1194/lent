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
    Foam::analyticalEllipse

SourceFiles
    analyticalEllipse.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Specialization of the analyticalSurface class for an ellipse.

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

#include "analyticalEllipse.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalEllipse, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalEllipse, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
vector analyticalEllipse::computeOneBySemiAxisSquare(const vector& semiAxes) const
{
    vector oneBySemiAxisSquare{};

    forAll(semiAxes, I)
    {
        oneBySemiAxisSquare[I] = 1.0/(semiAxes[I]*semiAxes[I]);
    }

    // This gives zero for the empty direction and is equivalent
    // to an infinite semi axis: this transforms the ellipsoid
    // into an ellipse
    oneBySemiAxisSquare = projector_&oneBySemiAxisSquare;

    return oneBySemiAxisSquare;
}

point analyticalEllipse::projectToReferencePlane(const point& aPoint) const
{
    return projector_&aPoint;
}

point analyticalEllipse::emptyComponent(const point& aPoint) const
{
    return ((Identity<scalar>{} - projector_) & aPoint);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalEllipse::analyticalEllipse(const dictionary& configDict)
:
    analyticalEllipsoid{configDict},
    projector_{}
{
    vector emptyDirection = configDict.lookup("emptyDirection");
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;

    vector centre{configDict.lookup("centre")};
    centre = projector_&centre;
    setCentre(centre);

    vector semiAxes{configDict.lookup("semiAxes")};
    auto oneBySemiAxisSquare = computeOneBySemiAxisSquare(semiAxes);
    setOneBySemiAxisSquare(oneBySemiAxisSquare);
}

analyticalEllipse::analyticalEllipse(const point& centre, const vector& semiAxes, const vector& emptyDirection)
:
   analyticalEllipsoid{centre, semiAxes},
   projector_{}
{
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;

    auto projectedCentre = projector_&centre;
    setCentre(projectedCentre);

    auto oneBySemiAxisSquare = computeOneBySemiAxisSquare(semiAxes);
    setOneBySemiAxisSquare(oneBySemiAxisSquare);
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalEllipse::distance(const point& trialPoint) const
{
    auto projection = projectToReferencePlane(trialPoint);

    return analyticalEllipsoid::distance(projection);
}

scalar analyticalEllipse::signedDistance(const point& trialPoint) const
{
    auto projection = projectToReferencePlane(trialPoint);

    return analyticalEllipsoid::signedDistance(projection);
}

point analyticalEllipse::normalProjectionToSurface(point& trialPoint) const
{
    auto projection = projectToReferencePlane(trialPoint);

    return analyticalEllipsoid::normalProjectionToSurface(projection)
            + emptyComponent(trialPoint);
}

vector analyticalEllipse::normalToPoint(const point& trialPoint) const
{
    auto projection = projectToReferencePlane(trialPoint);

    return analyticalEllipsoid::normalToPoint(projection);
}

point analyticalEllipse::intersection(const point& pointA,
                                     const point& pointB) const
{
    auto projectionA = projectToReferencePlane(pointA);
    auto projectionB = projectToReferencePlane(pointB);

    auto interSec = analyticalEllipsoid::intersection(projectionA, projectionB);

    // Sweet linearity: use the ratio of the projected points and their
    // intersection to compute a distance ratio. Then contruct the
    // intersection of the original points with this distance ratio.
    auto distanceBA = mag(projectionA - projectionB);
    auto distanceBIntersection = mag(interSec - projectionB);
    auto lambda = distanceBIntersection/distanceBA;

    return pointB + lambda*(pointA - pointB);
}

scalar analyticalEllipse::curvatureAt(const point& p) const
{
    auto projection = projectToReferencePlane(p);

    return analyticalEllipsoid::curvatureAt(projection);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
analyticalEllipse& analyticalEllipse::operator=(const analyticalEllipse& rhs)
{
    if (this != &rhs)
    {
        projector_ = rhs.projector_;
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
