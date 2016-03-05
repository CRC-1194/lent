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

    scalar projectedDistanceToOrigin = unitNormal_ & trialPoint;

    projected = projectedDistanceToOrigin / distanceOrigin_ * trialPoint;

    return projected;
}

vector analyticalPlane::normalToPoint(const point& trialPoint) const
{
    return unitNormal_;
}

//TODO: for intersection: use distance weighted blending of edge points
//      to obtain intersection with plane (analogue normalProjection)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
