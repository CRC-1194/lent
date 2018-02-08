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
    Foam::analyticalEllipsoid

SourceFiles
    analyticalEllipsoid.C

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

#include "analyticalEllipsoid.H"
#include "addToRunTimeSelectionTable.H"

#include <assert.h>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalEllipsoid, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalEllipsoid, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
scalar analyticalEllipsoid::levelSetValueOf(const point& aPoint) const
{
    // This function resembles
    // f(aPoint) = (x/a)^2 + (y/b)^2 + (z/c)^2 - 1
    scalar levelSetValue = -1.0;

    forAll(aPoint, I)
    {
        levelSetValue += pow(aPoint[I]/semiAxes_[I],2.0);
    }

    return levelSetValue; 
}

vector analyticalEllipsoid::levelSetGradientAt(const point& aPoint) const
{
    vector levelSetGradient{0,0,0};

    if (mag(aPoint - centre_) < SMALL)
    {
        FatalErrorInFunction
            << "Cannot compute level set gradient for a point which "
            << " coincides with the ellipsoid centre."
            << abort(FatalError);
    }

    forAll(levelSetGradient, I)
    {
        levelSetGradient[I] = 2.0/pow(semiAxes_[I],2.0)*aPoint[I];
    }

    return levelSetGradient;
}

point analyticalEllipsoid::moveToReferenceFrame(const point& aPoint) const
{
    return aPoint - centre_;
}

analyticalEllipsoid::parameterPair analyticalEllipsoid::intersectEllipsoidWithLine(const point& refPoint, const vector& path) const
{
    // Idea for intersection: inserting the parametric description of a line
    // composed by a refPoint and its path (direction) results in a quadratic
    // equation for the intersection parameter lambda.
    parameterPair lambdas{};

    scalar a = 0.0;
    scalar b = 0.0;
    scalar c = -1.0;

    forAll(refPoint, I)
    {
        a += std::pow(path[I]/semiAxes_[I],2.0);
        b += 2.0*refPoint[I]*path[I]/std::pow(semiAxes_[I],2.0);
        c += std::pow(refPoint[I]/semiAxes_[I],2.0);
    }

    // p-q formula...
    auto p = b/a;
    auto q = c/a;

    // Since the directional vector is normal to the ellipse
    // there must be two real valued solutions for the factor lambda
    // corresponding to the two intersections of the line with the
    // ellipse
    lambdas[0] = -0.5*p + std::sqrt(std::pow(0.5*p,2.0) - q);
    lambdas[1] = -0.5*p - std::sqrt(std::pow(0.5*p,2.0) - q);

    return lambdas;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalEllipsoid::analyticalEllipsoid(const dictionary& configDict)
:
    analyticalSurface{configDict},
    centre_{configDict.lookup("centre")},
    semiAxes_{configDict.lookup("semiAxes")}
{
}

analyticalEllipsoid::analyticalEllipsoid(const point& centre, const vector& semiAxes, const scalar& tolerance)
:
   analyticalSurface{},
   centre_{centre},
   semiAxes_{semiAxes}
{
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalEllipsoid::distance(const point& trialPoint) const
{
    return fabs(signedDistance(trialPoint));
}

scalar analyticalEllipsoid::signedDistance(const point& trialPoint) const
{
    auto copyPoint{trialPoint};
    auto projectedPoint = moveToReferenceFrame(trialPoint);
    auto surfaceProjection = normalProjectionToSurface(copyPoint);

    return sign(levelSetValueOf(projectedPoint))*mag(surfaceProjection - projectedPoint);
}

point analyticalEllipsoid::normalProjectionToSurface(point& trialPoint) const
{
    // Idea to find the normal projection on the ellipsoid: set up a line
    // using the refPoint and the level set gradient (which is normal to
    // the ellipsoid) as direction. Intersect this line with the level set
    // description of the ellipse which results in a quadratic equation to
    // be solved
    auto refPoint = moveToReferenceFrame(trialPoint);

    auto lGrad = levelSetGradientAt(refPoint);

    auto lambdas = intersectEllipsoidWithLine(refPoint, lGrad);

    // We want the intersection closer to refPoint, thus the lambda
    // with the smaller absolute value is the one we are looking for
    scalar lambda = 0.0;
    mag(lambdas[0]) < mag(lambdas[1]) ? lambda = lambdas[0] : lambda = lambdas[1];

    point ellipsePoint = refPoint + lambda*lGrad;

    assert(mag(levelSetValueOf(ellipsePoint)) < SMALL 
            && "Bug: projected point is not on ellipse");

    // Reverse movement to reference system
    return ellipsePoint + centre_; 
}

vector analyticalEllipsoid::normalToPoint(const point& trialPoint) const
{
    auto refPoint = moveToReferenceFrame(trialPoint);
    auto gradient = levelSetGradientAt(refPoint);
    auto levelSetValue = levelSetValueOf(refPoint);

    return sign(levelSetValue)*gradient/mag(gradient);
}

point analyticalEllipsoid::intersection(const point& pointA,
                                     const point& pointB) const
{
    // Move to reference frame
    auto refPointA = pointA - centre_;
    auto refPointB = pointB - centre_;

    auto lSetValueA = levelSetValueOf(refPointA);
    auto lSetValueB = levelSetValueOf(refPointB);

    // No intersection with elipsoid between A and B if their level set value
    // has the same sign
    // TODO: special case in which A and B form a tangent to the ellipsoid
    // is not covered yet. (TT)
    if (sign(lSetValueA) == sign(lSetValueB))
    {
        // Fatal Error
        FatalErrorInFunction
            << "Error: there is no intersection with ellipsoid between "
            << "point " << pointA << " and point " << pointB
            << abort(FatalError);
    }
    else if (lSetValueA == 0.0)
    {
        return pointA;
    }
    else if (lSetValueB == 0.0)
    {
        return pointB;
    }

    auto connection = refPointA - refPointB;
    auto basePoint = refPointB;

    // Ensure the direction points from inside the ellipsoid to the outside.
    // In this case solution of the quadratic equation will yield a positive
    // and a negative lambda with the positive one being the sought after (TT)
    if (lSetValueA < 0.0)
    {
        connection = refPointB -refPointA;
        basePoint = refPointA;
    }
    
    // Compute intersection using line approach
    auto lamdbas = intersectEllipsoidWithLine(refPointA, connection);
    
    // Decide which lambda is the sought after
    scalar lamdba = lamdbas[0];

    if (lamdbas[1] > 0.0)
    {
        lamdba = lamdbas[1];
    }
    
    return basePoint + lamdba*connection + centre_;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
analyticalEllipsoid& analyticalEllipsoid::operator=(const analyticalEllipsoid& rhs)
{
    if (this != &rhs)
    {
        centre_ = rhs.centre_;
        semiAxes_= rhs.semiAxes_;
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
