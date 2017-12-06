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

namespace Foam {
namespace FrontTracking {

    using constant::mathematical::pi;

    defineTypeNameAndDebug(analyticalEllipsoid, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalEllipsoid, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
scalar analyticalEllipsoid::devCosine(const scalar& longitude, const scalar& latitude, const point& p) const
{
    vector tLongitude{
                -semiAxes_.x()*sin(longitude)*sin(latitude),
                semiAxes_.y()*cos(longitude)*sin(latitude),
                0.0
            };

    vector tLatitude{
                semiAxes_.x()*cos(longitude)*cos(latitude),
                semiAxes_.y()*sin(longitude)*cos(latitude),
                -semiAxes_.z()*sin(latitude)
            };

    auto normal = tLatitude ^ tLongitude;
    auto connection = p - ellipsoidPoint(longitude, latitude);

    return mag(normal&connection) / (mag(normal)*mag(connection) + SMALL);
}

bool analyticalEllipsoid::converged(const scalar& longitude, const scalar& latitude, const point& p) const
{
    auto dev = devCosine(longitude, latitude, p);
    auto dist = distanceFromParameters(longitude, latitude, p);
    auto tanDev = dist*sqrt(1.0 - dev*dev);

    if (tanDev/semiAxes_.z() < tolerance_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

point analyticalEllipsoid::ellipsoidPoint(const scalar& longitude, const scalar& latitude) const
{
    return vector{
                    semiAxes_.x()*cos(longitude)*sin(latitude),
                    semiAxes_.y()*sin(longitude)*sin(latitude),
                    semiAxes_.z()*cos(latitude)
                 };
}

scalar analyticalEllipsoid::distanceFromParameters(const scalar& longitude, const scalar& latitude, const point& p) const
{
    return mag(ellipsoidPoint(longitude, latitude) - p);
}

scalar analyticalEllipsoid::curvature(const scalar& longitude, const scalar& latitude) const
{
    // Define shorter aliases to increase readability of curvature formula...
    const auto& a = semiAxes_.x();
    const auto& b = semiAxes_.y();
    const auto& c = semiAxes_.z();

    const auto& u = longitude;
    const auto& v = latitude;

    const auto numerator =
        a*b*c*(3*(a*a + b*b) + 2.0*c*c + (a*a + b*b - 2.0*c*c)*cos(2*v) - 2.0*(a*a - b*b)*cos(2*u)*sin(v)*sin(v));
    const auto denominator =
        8.0*pow((a*a*b*b*cos(v)*cos(v) + c*c*(b*b*cos(u)*cos(u) + a*a*sin(u)*sin(u))*sin(v)*sin(v)), 1.5);

    // Assumes front normal point outwards --> curvature is negative
    // The factor 2.0 arises since we actually need twice the mean curvature
    return -2.0*numerator/denominator;
}

analyticalEllipsoid::interval analyticalEllipsoid::findParameters(const point& p) const
{
    interval uInterval{0.0, pi/2.0};
    interval vInterval{0.0, pi/2.0};

    scalar minDistance = GREAT;
    scalar currentDistance = 0.0;

    scalar u;
    scalar v;
    scalar increment;

    interval uMinInterval{0.0, pi/2.0};
    interval vMinInterval{0.0, pi/2.0};

    for (label I = 0; I < 5; ++I)
    {
        increment = (uInterval[1] - uInterval[0])/30.0;
        u = uInterval[0] + 0.5*increment;

        while (u < uInterval[1])
        {
            v = vInterval[0] + 0.5*increment;

            while (v < vInterval[1])
            {
                currentDistance = distanceFromParameters(u, v, p);

                if (currentDistance < minDistance)
                {
                    minDistance = currentDistance;
                    uMinInterval = interval{u - increment, u + increment};
                    vMinInterval = interval{v - increment, v + increment};
                }

                v += increment;
            }

            u += increment;
        }

        uInterval = uMinInterval;
        vInterval = vMinInterval;

        if (converged((uInterval[0] + uInterval[1])/2.0, (vInterval[0] + vInterval[1])/2.0, p))
        {
            break;
        }
    }

    return interval{(uInterval[0] + uInterval[1])/2.0, (vInterval[0] + vInterval[1])/2.0};
}

point analyticalEllipsoid::projectToPositiveQuadrant(const point& p) const
{
    point projection{p};

    forAll(projection, I)
    {
        projection[I] = mag(projection[I]);
    }

    return projection;
}

void analyticalEllipsoid::undoQuadrantProjection(point& projection, const point& refPoint) const
{
    forAll(projection, I)
    {
        projection[I] *= sign(refPoint[I]);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalEllipsoid::analyticalEllipsoid(const dictionary& configDict)
:
    analyticalSurface{configDict},
    centre_{configDict.lookup("centre")},
    semiAxes_{configDict.lookup("semiAxes")},
    tolerance_{readScalar(configDict.lookup("tolerance"))}
{
}

analyticalEllipsoid::analyticalEllipsoid(const point& centre, const vector& semiAxes, const scalar& tolerance)
:
   analyticalSurface{},
   centre_{centre},
   semiAxes_{semiAxes},
   tolerance_{tolerance}
{
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar analyticalEllipsoid::distance(const point& trialPoint) const
{
    return fabs(signedDistance(trialPoint));
}

scalar analyticalEllipsoid::signedDistance(const point& trialPoint) const
{
    auto copy = trialPoint;
    auto closestPointOnSurface = normalProjectionToSurface(copy);
    vector connection = trialPoint - closestPointOnSurface;
    vector centreToSurface = closestPointOnSurface - centre_;

    return (mag(connection)*sign(connection&centreToSurface));
}

point analyticalEllipsoid::normalProjectionToSurface(point& trialPoint) const
{
    // For consistency, move point to reference ellipsoid with centre 
    // at (0 0 0)
    trialPoint -= centre_;
    auto projection = projectToPositiveQuadrant(trialPoint);
    auto parameters = findParameters(projection);
    auto pointOnEllipsoid = ellipsoidPoint(parameters[0], parameters[1]);
    undoQuadrantProjection(pointOnEllipsoid, trialPoint);

    return pointOnEllipsoid + centre_;
}

vector analyticalEllipsoid::normalToPoint(const point& trialPoint) const
{
    auto copy = trialPoint;
    auto closestPointOnSurface = normalProjectionToSurface(copy);
    vector connection = trialPoint - closestPointOnSurface;

    return connection/(mag(connection)+SMALL);
}

point analyticalEllipsoid::intersection(const point& pointA,
                                     const point& pointB) const
{
    notImplemented("analyticalEllipsoid::intgersection(...)");
    return point{0,0,0};
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
