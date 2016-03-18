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

#include "analyticalPlane.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalPlane, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalPlane, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalPlane::analyticalPlane(const dictionary& configDict)
:
    analyticalSurface(configDict)
{
    refPoint_ = configDict.lookup("referencePoint");
    unitNormal_ = configDict.lookup("normalVector");
    unitNormal_ = unitNormal_ / mag(unitNormal_);
    distanceOrigin_ = unitNormal_ & refPoint_;
}

analyticalPlane::analyticalPlane(const point& refPoint, const vector& normal)
:
    analyticalSurface()
{
    refPoint_ = refPoint;
    unitNormal_ = normal;
    unitNormal_ /= mag(unitNormal_);
    distanceOrigin_ = unitNormal_ & refPoint_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalPlane::distance(const point& trialPoint) const
{
    return fabs(signedDistance(trialPoint));
}

scalar analyticalPlane::signedDistance(const point& trialPoint) const
{
    return (trialPoint & unitNormal_) - distanceOrigin_;
}

point analyticalPlane::normalProjectionToSurface(point& trialPoint) const
{
    point projected(0.0, 0.0, 0.0);

    // Projection method fails for zero vector, thus provide 
    // alternative method for this case
    if (fabs(trialPoint & trialPoint) > SMALL)
    {
        scalar projectedDistanceToOrigin = unitNormal_ & trialPoint;

        projected = projectedDistanceToOrigin / distanceOrigin_ * trialPoint;
    }
    else
    {
        projected = distanceOrigin_ * unitNormal_;
    }

    return projected;
}

vector analyticalPlane::normalToPoint(const point& trialPoint) const
{
    return unitNormal_;
}

point analyticalPlane::intersection(const point& pointA, const point& pointB) const
{
        point intersect(0.0, 0.0, 0.0);

        scalar distanceRatio = distance(pointA) / (distance(pointA)
                                 + distance(pointB));

        // In essence, use (distance weighted) central differencing scheme to
        // find intersetion with plane
        intersect = distanceRatio*pointB + (1.0 - distanceRatio) * pointA;

        return intersect;
}


// * * * * * * * * * * * * * * Member Operators* * * * * * * * * * * * * * //
analyticalPlane& analyticalPlane::operator=(const analyticalPlane& plane) 
{
    if (this != &plane)
    {
        refPoint_ = plane.refPoint_;
        unitNormal_ = plane.unitNormal_;
    }

    return *this;
}


// * * * * * * * * * * * * * * Self Test * * * * * * * * * * * * * * * * * * //
void analyticalPlane::selfTest()
{
    Info << "\nStarting selftest of class analyticalPlane...\n" << endl;

    point testPointA(0.0, 0.0, 0.0);
    point testPointB(1000.0, 500.0, 1000.0);

    Info << "Trial point A: " << testPointA
         << "; trial point B: " << testPointB
         << endl;
    Info << "Signed distance: " << signedDistance(testPointA) << endl;
    Info << "Normal projection: " << normalProjectionToSurface(testPointA)
         << endl;
    Info << "Normal to point: " << normalToPoint(testPointA) << endl;
    Info << "Intersection: " << intersection(testPointA, testPointB) << endl;
    
    Info << "Finished tests\n" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
