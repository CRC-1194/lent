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

#include <assert.h>
#include <cmath>
#include <iomanip>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalEllipse, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalEllipse, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void analyticalEllipse::ensureValidHalfAxes()
{
    assert(semiAxes_.x() > 0.0 && semiAxes_.y() > 0.0 && semiAxes_.z() > 0.0
            && "Error: use semi axes larger than zero.");

    // For this class to work correctly the half axis value of the
    // empty direction can be anything but zero
    // In order to identify the empty axis, its value is set to the
    // non-sense value of -1 (TT)
    vector emptySemiAxis{-1, -1, -1};

    semiAxes_ = (projector_&semiAxes_)
                + ((Identity<scalar>{} - projector_)&emptySemiAxis);
}

void analyticalEllipse::ensureValidCentre()
{
    // The empty direction component of the centre has to be zero
    // for this class to work properly
    centre_ = projector_&centre_;
}

label analyticalEllipse::majorSemiAxisIndex() const
{
    // FIXME: for some reason the OpenFOAM version of 'max()' does
    // not work here (TT)
    label indexMax = 0;
        
    forAll(semiAxes_, I)
    {
        if (semiAxes_[I] > semiAxes_[indexMax])
        {
            indexMax = I;
        }
    }

    return indexMax;
}

label analyticalEllipse::minorSemiAxisIndex() const
{
    auto majorIndex = majorSemiAxisIndex();

    forAll(semiAxes_, I)
    {
        if (I != majorIndex && semiAxes_[I] > 0.0)
        {
            return I;
        }
    }

    // This should never be reached (TT)
    return 0;
}

scalar analyticalEllipse::levelSetValueOf(const point& aPoint) const
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

vector analyticalEllipse::levelSetGradientAt(const point& aPoint) const
{
    vector levelSetGradient{0,0,0};

    if (mag(aPoint - centre_) < SMALL)
    {
        FatalErrorInFunction
            << "Cannot compute level set gradient for a point which "
            << " coincides with the ellipse centre."
            << abort(FatalError);
    }

    forAll(levelSetGradient, I)
    {
        levelSetGradient[I] = 2.0/pow(semiAxes_[I],2.0)*aPoint[I];
    }

    return levelSetGradient;
}

point analyticalEllipse::moveToReferenceFrame(const point& aPoint) const
{
    return (projector_&aPoint) - centre_;
}

vector analyticalEllipse::emptyComponent(const point& aPoint) const
{
    return ((Identity<scalar>{} - projector_) & aPoint);
}

analyticalEllipse::parameterPair analyticalEllipse::intersectEllipseWithLine(const point& refPoint, const vector& path) const
{
    parameterPair lambdas{};

    // Idea for intersection: inserting the parametric description of a line
    // composed by a refPoint and its path (direction) results in a quadratic
    // equation for the intersection parameter lambda.
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

    // It is up to the user of this function to ensure
    // that an intersection exists. Thus a negative discriminant
    // is considered an error.
    assert((0.25*p*p - q) >= 0.0 &&
            "Error in analyticalEllipse::intersectEllipseWithLine: no intersection with ellipse found.");

    lambdas[0] = -0.5*p + std::sqrt(std::pow(0.5*p,2.0) - q);
    lambdas[1] = -0.5*p - std::sqrt(std::pow(0.5*p,2.0) - q);

    return lambdas;
}

scalar analyticalEllipse::ellipseCurvature(const point& aPoint) const
{
    scalar curvature = 0.0;

    vector oneBySemiAxisSqr{0,0,0};

    forAll(oneBySemiAxisSqr, I)
    {
        oneBySemiAxisSqr[I] = 1.0/(semiAxes_[I]*semiAxes_[I] + SMALL);
    }

    // Remove empty component
    oneBySemiAxisSqr = projector_&oneBySemiAxisSqr;

    auto g = levelSetGradientAt(aPoint);
    auto magGrad = mag(g);
    
    curvature = oneBySemiAxisSqr.x()*(g.y()*g.y() + g.z()*g.z())
              + oneBySemiAxisSqr.y()*(g.x()*g.x() + g.z()*g.z())
              + oneBySemiAxisSqr.z()*(g.x()*g.x() + g.y()*g.y());

    return -2.0/(magGrad*magGrad*magGrad)*curvature;
}

scalar analyticalEllipse::bisection(const std::function< scalar(scalar)>& rootFunction, analyticalEllipse::parameterPair interval) const
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
analyticalEllipse::analyticalEllipse(const dictionary& configDict)
:
    analyticalSurface{configDict},
    centre_{configDict.lookup("centre")},
    semiAxes_{configDict.lookup("semiAxes")}
{
    vector emptyDirection = configDict.lookup("emptyDirection");
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;
    ensureValidHalfAxes();
    ensureValidCentre();
}

analyticalEllipse::analyticalEllipse(const point& centre, const vector& semiAxes, const vector& emptyDirection)
:
   analyticalSurface{},
   centre_{centre},
   semiAxes_{semiAxes}
{
    projector_ = Identity<scalar>{} - emptyDirection*emptyDirection;
    ensureValidHalfAxes();
    ensureValidCentre();
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalEllipse::distance(const point& trialPoint) const
{
    return fabs(signedDistance(trialPoint));
}

scalar analyticalEllipse::signedDistance(const point& trialPoint) const
{
    auto copyPoint{trialPoint};
    auto projectedPoint = moveToReferenceFrame(trialPoint);
    auto surfaceProjection = normalProjectionToSurface(copyPoint);

    return sign(levelSetValueOf(projectedPoint))*mag(surfaceProjection - projectedPoint);
}

point analyticalEllipse::normalProjectionToSurface(point& trialPoint) const
{
    point pointOnSurface{0,0,0};

    // Move trial point to first quadrant of reference frame
    auto refPoint = moveToReferenceFrame(trialPoint);
    auto p = refPoint;
    forAll(p, I)
    {
        p[I] = fabs(p[I]);
    }

    vector a{semiAxes_};
    auto minDistanceParameter = [&a,&p](scalar t)
            {
                auto term1 = a.x()*p.x()/(a.x()*a.x() + t);
                auto term2 = a.y()*p.y()/(a.y()*a.y() + t);
                auto term3 = a.z()*p.z()/(a.z()*a.z() + t);

                return term1*term1 + term2*term2 + term3*term3 - 1.0;
            };

    // The minor and the major semi axis are rquired to compute the
    // appropriate parameter range (TT)
    parameterPair interval{};
    auto minI = minorSemiAxisIndex();
    auto maxI = majorSemiAxisIndex();
    interval[0] = -1.0*a[minI]*a[minI] + a[minI]*p[minI];
    interval[1] = -1.0*a[minI]*a[minI]
                 + sqrt(a[maxI]*a[maxI]*p[maxI]*p[maxI]
                        + a[minI]*a[minI]*p[minI]*p[minI]);

    auto lambdaMin = bisection(minDistanceParameter, interval);
    
    pointOnSurface.x() = a.x()*a.x()*p.x() / (a.x()*a.x() + lambdaMin);
    pointOnSurface.y() = a.y()*a.y()*p.y() / (a.y()*a.y() + lambdaMin);
    pointOnSurface.z() = a.z()*a.z()*p.z() / (a.z()*a.z() + lambdaMin);

    // Move point on surface to correct quadrant
    forAll(pointOnSurface, I)
    {
        pointOnSurface[I] *= sign(refPoint[I]);
    }

    return pointOnSurface + centre_ + emptyComponent(trialPoint);
}

vector analyticalEllipse::normalToPoint(const point& trialPoint) const
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

point analyticalEllipse::intersection(const point& pointA,
                                     const point& pointB) const
{
    auto refPointA = moveToReferenceFrame(pointA);
    auto refPointB = moveToReferenceFrame(pointB);

    auto lSetValueA = levelSetValueOf(refPointA);
    auto lSetValueB = levelSetValueOf(refPointB);

    // No intersection with elipsoid between A and B if their level set value
    // has the same sign
    // TODO: special case in which A and B form a tangent to the ellipse
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
    auto emptyComponent = (Identity<scalar>{} - projector_)&pointB;

    // Ensure the direction points from inside the ellipsoid to the outside.
    // In this case solution of the quadratic equation will yield a positive
    // and a negative lambda with the positive one being the sought after (TT)
    if (lSetValueA < 0.0)
    {
        connection = refPointB -refPointA;
        basePoint = refPointA;
        emptyComponent = (Identity<scalar>{} - projector_)&pointA;
    }
    
    // Compute intersection using line approach
    auto lamdbas = intersectEllipseWithLine(refPointA, connection);
    
    // Decide which lambda is the sought after
    scalar lamdba = lamdbas[0];

    if (lamdbas[1] > 0.0)
    {
        lamdba = lamdbas[1];
    }
    
    return basePoint + lamdba*connection + centre_ + emptyComponent;
}

scalar analyticalEllipse::curvatureAt(const point& aPoint) const
{
    auto refPoint = moveToReferenceFrame(aPoint);

    if (mag(levelSetValueOf(refPoint)) > SMALL)
    {
        // Remember: the normal projection moves the point back to the
        // original coordinate system
        point copyPoint{aPoint};
        refPoint = normalProjectionToSurface(copyPoint);
        refPoint = moveToReferenceFrame(refPoint);
    }

    return ellipseCurvature(refPoint); 
}

void analyticalEllipse::centre(const point& newCentre)
{
    centre_ = newCentre;
    ensureValidCentre(); 
}

void analyticalEllipse::semiAxes(const vector& newSemiAxes)
{
    semiAxes_ = newSemiAxes;
    ensureValidHalfAxes();
}

void analyticalEllipse::writeParameters(const word fileName) const
{
    // TODO: open in append mode (TT)
    OFstream outputFile(fileName);

    outputFile.stdStream() << std::setprecision(15);

    outputFile << "-------------------------------\n"
               << "type " << this->type() << '\n'
               << "centre " << centre_ << '\n'
               << "semiAxes " << semiAxes_ << '\n';
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
analyticalEllipse& analyticalEllipse::operator=(const analyticalEllipse& rhs)
{
    if (this != &rhs)
    {
        centre_ = rhs.centre_;
        semiAxes_ = rhs.semiAxes_;
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
