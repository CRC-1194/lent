/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "analyticalSphere.H"
#include "addToRunTimeSelectionTable.H"

#include <cmath>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalSphere, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalSphere, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalSphere::analyticalSphere(const dictionary& configDict)
:
    analyticalSurface(configDict)
{
    centre_ = configDict.lookupOrDefault<vector>("centre", centre_);
    radius_ = configDict.lookupOrDefault<scalar>("radius", radius_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalSphere::distance(const point& trialPoint) const
{
    return fabs(signedDistance(trialPoint));
}

scalar analyticalSphere::signedDistance(const point& trialPoint) const
{
    return mag(trialPoint - centre_) - radius_;
}

point analyticalSphere::normalProjectionToSurface(point& trialPoint) const
{
    point projected(0.0, 0.0, 0.0);
    vector normalizedDirection = (trialPoint - centre_)
                                    / mag(trialPoint - centre_);
    return radius_*normalizedDirection + centre_;
}

vector analyticalSphere::normalToPoint(const point& trialPoint) const
{
    vector normal(trialPoint - centre_);
    return normal/mag(normal);
}

point analyticalSphere::intersection(const point& pointA,
                                     const point& pointB) const
{
    // Basic idea for intersection:
    // 1) reduce to planar problem by intersecting the sphere with the plane
    //      (centre_, pointA, pointB)
    //      --> results in circle with same radius as sphere
    // 2) exploit fact that the intersection will always have a distance 
    //      to the centre equivalent to the radius 
    // 3) use trigenometry to determine the distance ratio in which the segement
    //      pointA--pointB is intersected
    // 4) compute intersection using inverse distance weighted interpolation
    point intersect(0.0, 0.0, 0.0);

    // TODO: determine which point lies inside the sphere and ensure that 
    // it is used inthe following
    scalar distanceCA = mag(pointA - centre_);
    scalar distanceAB = mag(pointA - pointB);

    // Angle enclosed by centre_ --> pointA and centre_ --> intersection
    scalar alpha = std::acos(distanceCA / radius_);

    scalar distanceRatio = std::sin(alpha)*radius_ / distanceAB;

    intersect = distanceRatio*pointB + (1.0 - distanceRatio)*pointA;

    return intersect;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
