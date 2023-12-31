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
    Foam::FrontTracking:analyticalSphere

SourceFiles
    analyticalSphere.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for a sphere.

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

#include "analyticalSphere.H"
#include "addToRunTimeSelectionTable.H"

#include <iomanip>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalSphere, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalSphere, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalSphere::analyticalSphere(const dictionary& configDict)
:
    analyticalSurface{configDict}
{
    centre_ = configDict.get<point>("centre");
    radius_ = configDict.get<scalar>("radius");
    name_ = configDict.lookupOrDefault<word>("name", analyticalSphere::typeName);
}

analyticalSphere::analyticalSphere(
    const point& centre, 
    const scalar radius,
    const word name 
)
:
   analyticalSurface{},
   centre_{centre}, 
   radius_{radius}, 
   name_{name}
{} 


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
    // Note: projection does not work if trialPoint coincides with
    //       centre_
    point projected(0.0, 0.0, 0.0);

    // Ensure sufficient distance between trialPoint and centre_
    if (magSqr(trialPoint - centre_) > SMALL)
    {
        vector normalizedDirection = (trialPoint - centre_)
                                      / mag(trialPoint - centre_);
        return radius_*normalizedDirection + centre_;
    }
    else
    {
        Info << "Warning: point " << trialPoint << "coincides with "
             << "centre.\n"
             << "Cannot project, using centre instead" << endl;
        return centre_;
    }
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

    // Determine inner and outer point to ensure correctness of 
    // intersection method
    point innerPoint = pointA;
    point outerPoint = pointB;

    if (signedDistance(pointA) > 0.0)
    {
        innerPoint = pointB;
        outerPoint = pointA;
    } 

    // Catch exceptional case in which the inner point coincides with
    // the centre_ and aforementioned algorithm fails
    if (mag(innerPoint - centre_) > SMALL)
    {
        scalar distanceCI = mag(centre_ - innerPoint);
        scalar distanceInOut = mag(innerPoint - outerPoint);
        scalar cosAlpha = (centre_ - innerPoint) & (outerPoint - innerPoint)
                            / (distanceCI * distanceInOut);

        // Exploit cosine rule for an arbitrary triangle and solve
        // for the sought-after distance between the innerPoint and the
        // intersection using the pq-formula. Note, that the factor 0.5
        // cancels out. Only addition of the root term yields geometrically
        // valid results.
        scalar p = -distanceCI * cosAlpha;
        scalar q = sqr(distanceCI) - sqr(radius_);
        scalar distanceIintersect = -p + sqrt(sqr(p) - q);
        scalar distanceRatio = distanceIintersect / distanceInOut;

        intersect = distanceRatio * outerPoint +
                        (1.0 - distanceRatio) * innerPoint;
    }
    else
    {
        intersect = centre_ + radius_ * (outerPoint - centre_)
                        / mag(outerPoint - centre_);
    }

    return intersect;
}

void analyticalSphere::centre(const point& newCentre)
{
    centre_ = newCentre;
}

void analyticalSphere::radius(const scalar newRadius)
{
    radius_ = newRadius;
}

void analyticalSphere::writeParameters(const word fileName) const
{
    auto outputFile = outputStream(fileName);

    outputFile.stdStream() << std::setprecision(15);

    outputFile << "-------------------------------\n"
               << "type " << this->type() << '\n'
               << "centre " << centre_ << '\n'
               << "radius " << radius_ << '\n';
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
analyticalSphere& analyticalSphere::operator=(const analyticalSphere& rhs)
{
    if (this != &rhs)
    {
        centre_ = rhs.centre_;
        radius_ = rhs.radius_;
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
