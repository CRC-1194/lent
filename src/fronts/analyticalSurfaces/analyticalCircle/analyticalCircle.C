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
    Foam::analyticalCircle

SourceFiles
    analyticalCircle.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Specialization of the analyticalSurface class for a sphere.

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

#include "analyticalCircle.H"
#include "addToRunTimeSelectionTable.H"

#include <cassert>
#include <iomanip>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalCircle, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalCircle, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalCircle::analyticalCircle(const dictionary& configDict)
:
    analyticalSurface(configDict)
{
    centre_ = configDict.lookupOrDefault<vector>("centre", centre_);
    radius_ = configDict.lookupOrDefault<scalar>("radius", radius_);
    vector emptyDirection = configDict.lookup("emptyDirection");
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;

    centre_ = projector_&centre_;
}

analyticalCircle::analyticalCircle(const point& centre, const scalar radius, const vector& emptyDirection)
:
   analyticalSurface{},
   centre_{centre},
   radius_{radius}
{
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;
    centre_ = projector_&centre_;
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalCircle::distance(const point& trialPoint) const
{
    return fabs(signedDistance(trialPoint));
}

scalar analyticalCircle::signedDistance(const point& trialPoint) const
{
    return mag((projector_&trialPoint) - centre_) - radius_;
}

point analyticalCircle::normalProjectionToSurface(point& trialPoint) const
{
    // Note: projection does not work if trialPoint coincides with
    //       centre_
    point projected(0.0, 0.0, 0.0);

    vector centreToPoint = projector_&(trialPoint - centre_);

    // Ensure sufficient distance between trialPoint and centre_
    if (magSqr(centreToPoint) > SMALL)
    {
        vector normalizedDirection = centreToPoint / mag(centreToPoint);
        return radius_*normalizedDirection + centre_
                    + ((Identity<scalar>{} - projector_)&trialPoint);
    }
    else
    {
        Info << "Warning: point " << trialPoint << "coincides with "
             << "centre.\n"
             << "Cannot project, using centre instead" << endl;
        return centre_ + ((Identity<scalar>{} - projector_)&trialPoint);
    }
}

vector analyticalCircle::normalToPoint(const point& trialPoint) const
{
    vector normal(projector_&(trialPoint - centre_));
    return normal/mag(normal);
}

point analyticalCircle::intersection(const point& pointA,
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
    point innerPoint = projector_&pointA;
    point outerPoint = projector_&pointB;

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

    // Add component in empty direction by linear interpolation
    auto distA = distance(pointA);
    auto distB = distance(pointB);
    auto distAB = distA + distB;

    vector emptyComponent = (Identity<scalar>{} - projector_)&
                    (distB/distAB*pointA + (1.0 - distB/distAB)*pointB);

    return intersect + emptyComponent;
}

void analyticalCircle::centre(const vector newCentre)
{
    centre_ = projector_&newCentre;
}

void analyticalCircle::radius(const scalar newRadius)
{
    assert(newRadius > 0.0 && "Error: radius of analyticalCircle must be greater than zero.");
    radius_ = newRadius;
}


void analyticalCircle::writeParameters(const word fileName) const
{
    auto outputFile = outputStream(fileName);

    outputFile.stdStream() << std::setprecision(15);

    outputFile << "-------------------------------\n"
               << "type " << this->type() << '\n'
               << "centre " << centre_ << '\n'
               << "radius " << radius_ << '\n';
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
analyticalCircle& analyticalCircle::operator=(const analyticalCircle& rhs)
{
    if (this != &rhs)
    {
        centre_ = rhs.centre_;
        radius_ = rhs.radius_;
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
