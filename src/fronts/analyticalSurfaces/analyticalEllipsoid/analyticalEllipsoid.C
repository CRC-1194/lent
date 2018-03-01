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
    Specialization of the analyticalSurface class for an ellipsoid.

    See https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
    for information on how to find the minimal distance between a point
    and the ellipsoid.

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

#include <cmath>
#include <assert.h>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalEllipsoid, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalEllipsoid, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
scalar analyticalEllipsoid::levelSetValueOf(const point& aPoint) const
{
    // This function resembles
    // f(aPoint) = qx^2 + ry^2 + sz^2 - 1, where q, r, s are one by semi-axis
    // squared
    scalar levelSetValue = -1.0;

    forAll(aPoint, I)
    {
        levelSetValue += oneBySemiAxisSqr_[I]*aPoint[I]*aPoint[I];
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
        levelSetGradient[I] = 2.0*oneBySemiAxisSqr_[I]*aPoint[I];
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

    scalar quadCoeff = 0.0;
    scalar linCoeff = 0.0;
    scalar absCoeff = -1.0; // Already include the "-1" from the level set equation

    forAll(refPoint, I)
    {
        quadCoeff += oneBySemiAxisSqr_[I]*path[I]*path[I];
        linCoeff += 2.0*oneBySemiAxisSqr_[I]*refPoint[I]*path[I];
        absCoeff += oneBySemiAxisSqr_[I]*refPoint[I]*refPoint[I];
    }

    // p-q formula...
    auto p = linCoeff/quadCoeff;
    auto q = absCoeff/quadCoeff;

    // It is the responsibility of the user of this class to request
    // existing intersections only. Thus, a negative discriminant will
    // terminate the program
    assert((0.25*linCoeff*linCoeff - absCoeff) >= 0.0 &&
            "analyticalEllipsoid::intersectEllipsoidWithLine: no real-valued intersection found.");

    lambdas[0] = -0.5*p + std::sqrt(std::pow(0.5*p,2.0) - q);
    lambdas[1] = -0.5*p - std::sqrt(std::pow(0.5*p,2.0) - q);

    return lambdas;
}

scalar analyticalEllipsoid::ellipsoidCurvature(const point& p) const
{
    // See the notes on the level set equation used for fitting for a
    // derivation of the ellipsoid curvature from a level set representation.
    // It boils down to write the level set equation of an ellipsoid
    // phi = qx^2 + ry^2 + sz^2 - 1 = 0 and apply
    // kappa = -div(grad(phi)/mag(grad(phi))).
    auto g = levelSetGradientAt(p);
    auto D = mag(g);

    scalar kappa = (g.y()*g.y() + g.z()*g.z())*oneBySemiAxisSqr_.x()
                 + (g.x()*g.x() + g.z()*g.z())*oneBySemiAxisSqr_.y()
                 + (g.x()*g.x() + g.y()*g.y())*oneBySemiAxisSqr_.z();

    return -2.0*kappa/(D*D*D + SMALL);
}

label analyticalEllipsoid::minorSemiAxisIndex() const
{
    label minIndex = 0;
    scalar minAxis = 0.0;

    forAll(oneBySemiAxisSqr_, I)
    {
        if (oneBySemiAxisSqr_[I] > minAxis)
        {
            minAxis = oneBySemiAxisSqr_[I];
            minIndex = I;
        }
    }

    return minIndex;
}

scalar analyticalEllipsoid::bisection(const std::function< scalar(scalar)>& rootFunction, analyticalEllipsoid::parameterPair interval) const
{
    auto midPoint = 0.5*(interval[0] + interval[1]);
    scalar valueAtMidPoint = rootFunction(midPoint);

    // TODO: set hard limit for iteration count? (TT)
    while(mag(valueAtMidPoint) > 1.0e-13 && mag(midPoint - interval[1]) > 1.0e-13)
    {
        if (rootFunction(interval[0])*valueAtMidPoint > 0.0)
        {
            interval[0] = midPoint;
        }
        else
        {
            interval[1] = midPoint;
        }

        midPoint = 0.5*(interval[0] + interval[1]);
        valueAtMidPoint = rootFunction(midPoint);
    }

    return midPoint;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalEllipsoid::analyticalEllipsoid(const dictionary& configDict)
:
    analyticalSurface{configDict},
    centre_{configDict.lookup("centre")},
    oneBySemiAxisSqr_{}
{
    vector semiAxes{configDict.lookup("semiAxes")};

    forAll(semiAxes, I)
    {
        oneBySemiAxisSqr_[I] = 1.0/(semiAxes[I]*semiAxes[I]);
    }
}

analyticalEllipsoid::analyticalEllipsoid(const point& centre, const vector& semiAxes)
:
   analyticalSurface{},
   centre_{centre},
   oneBySemiAxisSqr_{}
{
    forAll(semiAxes, I)
    {
        oneBySemiAxisSqr_[I] = 1.0/(semiAxes[I]*semiAxes[I]);
    }
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
    point pointOnSurface{0,0,0};

    // Move trial point to first quadrant of reference frame
    auto refPoint = moveToReferenceFrame(trialPoint);
    auto p = refPoint;
    forAll(p, I)
    {
        p[I] = fabs(p[I]);
    }

    vector q{oneBySemiAxisSqr_};
    auto minDistanceParameter = [&q,&p](scalar t)
            {
                auto term1 = p.x()/(1.0 + q.x()*t);
                auto term2 = p.y()/(1.0 + q.y()*t);
                auto term3 = p.z()/(1.0 + q.z()*t);

                return q.x()*term1*term1 + q.y()*term2*term2
                        + q.z()*term3*term3 - 1.0;
            };

    // The minor and the major semi axis are rquired to compute the
    // appropriate parameter range (TT)
    parameterPair interval{};
    auto minI = minorSemiAxisIndex();
    interval[0] = -1.0/q[minI] + p[minI]/sqrt(q[minI]);
    interval[1] = -1.0/q[minI]
                 + sqrt(p.x()*p.x()/q.x() + p.y()*p.y()/q.y()
                         + p.z()*p.z()/q.z());

    auto lambdaMin = bisection(minDistanceParameter, interval);
    
    pointOnSurface.x() = p.x() / (1.0 + q.x()*lambdaMin);
    pointOnSurface.y() = p.y() / (1.0 + q.y()*lambdaMin);
    pointOnSurface.z() = p.z() / (1.0 + q.z()*lambdaMin);

    // Move point on surface to correct quadrant
    forAll(pointOnSurface, I)
    {
        pointOnSurface[I] *= sign(refPoint[I]);
    }

    return pointOnSurface + centre_;
}

vector analyticalEllipsoid::normalToPoint(const point& trialPoint) const
{
    // NOTE: using the level set gradient rather than the connection
    // between the refPoint and pointOnSurface is more stable for points
    // close to the interface (TT)
    point copyPoint{trialPoint};
    auto pointOnSurface = normalProjectionToSurface(copyPoint);
    pointOnSurface = moveToReferenceFrame(pointOnSurface);
    auto gradient = levelSetGradientAt(pointOnSurface);

    return gradient/mag(gradient);
}

point analyticalEllipsoid::intersection(const point& pointA,
                                     const point& pointB) const
{
    // NOTE: this function is intended to find a unique intersection
    // between A and B. The case of A and B spanning the ellipsoid
    // is not covered yet. (TT)
    
    // Move to reference frame
    auto refPointA = moveToReferenceFrame(pointA);
    auto refPointB = moveToReferenceFrame(pointB);

    auto lSetValueA = levelSetValueOf(refPointA);
    auto lSetValueB = levelSetValueOf(refPointB);

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

scalar analyticalEllipsoid::curvatureAt(const point& p) const
{
    point pCopy{p};
    auto pOnEllipsoid = normalProjectionToSurface(pCopy);
    pOnEllipsoid = moveToReferenceFrame(pOnEllipsoid);

    return ellipsoidCurvature(pOnEllipsoid);
}

vector analyticalEllipsoid::semiAxes() const
{
    vector semiAxes{};

    forAll(semiAxes, I)
    {
        semiAxes[I] = sqrt(1.0/(oneBySemiAxisSqr_[I] + SMALL));
    }

    return semiAxes;
}

// * * * * * * * * * * * * * * Member Operators    * * * * * * * * * * * * * * //
analyticalEllipsoid& analyticalEllipsoid::operator=(const analyticalEllipsoid& rhs)
{
    if (this != &rhs)
    {
        centre_ = rhs.centre_;
        oneBySemiAxisSqr_= rhs.oneBySemiAxisSqr_;
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
