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
    Foam::FrontTracking:analyticalEllipsoid

SourceFiles
    analyticalEllipsoid.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for an ellipsoid.

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

#include "analyticalEllipsoid.H"
#include "addToRunTimeSelectionTable.H"

#include <assert.h>
#include <cmath>
#include <iomanip>

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

    for (direction I = 0; I != 3; ++I)
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

    for (direction I = 0; I != 3; ++I)
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

    for (direction I = 0; I != 3; ++I)
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

direction analyticalEllipsoid::minorSemiAxisIndex() const
{
    direction minIndex = 0;
    scalar minAxis = 0.0;

    for (direction I = 0; I != 3; ++I)
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

scalar analyticalEllipsoid::signedDistanceRefSytem(const point P) const
{
    auto surfaceProjectionP = normalProjectionToSurfaceRefSystem(P);

    return sign(levelSetValueOf(P))*mag(surfaceProjectionP - P);
}

point analyticalEllipsoid::normalProjectionToSurfaceRefSystem(point P) const
{
    point pointOnSurface{0,0,0};
    const double epsilon = std::numeric_limits<double>::epsilon();

    // Move point to first quadrant of reference frame
    auto p = P;
    for (direction I = 0; I != 3; ++I)
    {
        p[I] = fabs(p[I]);
    }

    vector q{oneBySemiAxisSqr_};

    auto minDistanceParameter = [&q,&p,&epsilon](scalar t)
            {
                auto term1 = p.x()/(1.0 + q.x()*t + epsilon);
                auto term2 = p.y()/(1.0 + q.y()*t + epsilon);
                auto term3 = p.z()/(1.0 + q.z()*t + epsilon);

                return q.x()*term1*term1 + q.y()*term2*term2
                        + q.z()*term3*term3 - 1.0;
            };

    // The minor and the major semi axis are required to compute the
    // appropriate parameter range (TT)
    parameterPair interval{};
    auto minI = minorSemiAxisIndex();
    interval[0] = -1.0/q[minI] + p[minI]/sqrt(q[minI]);
    interval[1] = -1.0/q[minI]
                 + sqrt(p.x()*p.x()/q.x() + p.y()*p.y()/q.y()
                         + p.z()*p.z()/q.z());

    auto lambdaMin = bisection(minDistanceParameter, interval);
    
    pointOnSurface.x() = p.x() / (1.0 + q.x()*lambdaMin + epsilon);
    pointOnSurface.y() = p.y() / (1.0 + q.y()*lambdaMin + epsilon);
    pointOnSurface.z() = p.z() / (1.0 + q.z()*lambdaMin + epsilon);

    // Move point on surface to correct quadrant
    for (direction I = 0; I != 3; ++I)
    {
        pointOnSurface[I] *= sign(P[I]);
    }

    return pointOnSurface;
}

vector analyticalEllipsoid::normalToPointRefSystem(const point P) const
{
    // NOTE: using the level set gradient rather than the connection
    // between the refPoint and pointOnSurface is more stable for points
    // close to the interface (TT)
    auto pointOnSurface = normalProjectionToSurfaceRefSystem(P);
    auto gradient = levelSetGradientAt(pointOnSurface);

    return gradient/mag(gradient);
}

point analyticalEllipsoid::intersectionRefSystem(const point pointA, const point pointB) const
{
    // NOTE: this function is intended to find a unique intersection
    // between A and B. The case of A and B spanning the ellipsoid
    // is not covered yet. (TT)
    auto lSetValueA = levelSetValueOf(pointA);
    auto lSetValueB = levelSetValueOf(pointB);

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

    auto connection = pointA - pointB;
    auto basePoint = pointB;

    // Ensure the direction points from inside the ellipsoid to the outside.
    // In this case solution of the quadratic equation will yield a positive
    // and a negative lambda with the positive one being the sought after (TT)
    if (lSetValueA < 0.0)
    {
        connection = pointB -pointA;
        basePoint = pointA;
    }
    
    // Compute intersection using line approach
    auto lamdbas = intersectEllipsoidWithLine(pointA, connection);
    
    // Decide which lambda is the sought after
    scalar lamdba = lamdbas[0];

    if (lamdbas[1] > 0.0)
    {
        lamdba = lamdbas[1];
    }
    
    return basePoint + lamdba*connection;
}

scalar analyticalEllipsoid::curvatureAtRefSystem(const point P) const
{
    auto pOnEllipsoid = normalProjectionToSurfaceRefSystem(P);

    return ellipsoidCurvature(pOnEllipsoid);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

analyticalEllipsoid::analyticalEllipsoid(const dictionary& configDict)
:
    analyticalSurface{configDict},
    centre_{configDict.get<point>("centre")},
    oneBySemiAxisSqr_{}
{
    vector semiAxes{configDict.get<vector>("semiAxes")};

    for(direction I = 0; I != 3; ++I)
    {
        oneBySemiAxisSqr_[I] = 1.0/(semiAxes[I]*semiAxes[I]);
    }

    name_ = configDict.lookupOrDefault<word>("name", analyticalEllipsoid::typeName);
}

analyticalEllipsoid::analyticalEllipsoid(
    const point& centre, 
    const vector& semiAxes, 
    const word name
)
:
   analyticalSurface{},
   centre_{centre},
   oneBySemiAxisSqr_{},
   name_(name)
{
    for (direction I = 0; I != 3; ++I)
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
    return signedDistanceRefSytem(moveToReferenceFrame(trialPoint));
}

point analyticalEllipsoid::normalProjectionToSurface(point& trialPoint) const
{
    return centre_ + normalProjectionToSurfaceRefSystem(moveToReferenceFrame(trialPoint));
}

vector analyticalEllipsoid::normalToPoint(const point& trialPoint) const
{
    return normalToPointRefSystem(moveToReferenceFrame(trialPoint));
}

point analyticalEllipsoid::intersection(const point& pointA,
                                     const point& pointB) const
{
    return centre_ + intersectionRefSystem(moveToReferenceFrame(pointA),
                                            moveToReferenceFrame(pointB));
}

scalar analyticalEllipsoid::curvatureAt(const point& p) const
{
    return curvatureAtRefSystem(moveToReferenceFrame(p));
}

point analyticalEllipsoid::centre() const
{
    return centre_;
}

vector analyticalEllipsoid::semiAxes() const
{
    vector semiAxes{};

    for (direction I = 0; I != 3; ++I)
    {
        semiAxes[I] = sqrt(1.0/(oneBySemiAxisSqr_[I] + SMALL));
    }

    return semiAxes;
}

void analyticalEllipsoid::oneBySemiAxisSquare(const vector& oneBySemiAxisSquare)
{
    oneBySemiAxisSqr_ = oneBySemiAxisSquare;
}

void analyticalEllipsoid::centre(const vector& newCentre)
{
    centre_ = newCentre;
}

void analyticalEllipsoid::semiAxes(const vector& newSemiAxes)
{
    for (direction I = 0; I != 3; ++I)
    {
        oneBySemiAxisSqr_[I] = 1.0/(newSemiAxes[I]*newSemiAxes[I]);
    }
}

void analyticalEllipsoid::writeParameters(const word fileName) const
{
    auto outputFile = outputStream(fileName);

    outputFile.stdStream() << std::setprecision(15);

    outputFile << "-------------------------------\n"
               << "type " << this->type() << '\n'
               << "centre " << centre_ << '\n'
               << "semiAxes " << semiAxes() << '\n';
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
