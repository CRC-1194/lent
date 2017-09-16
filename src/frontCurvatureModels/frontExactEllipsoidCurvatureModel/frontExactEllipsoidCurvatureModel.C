/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "frontExactEllipsoidCurvatureModel.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    using constant::mathematical::pi;

    defineTypeNameAndDebug(frontExactEllipsoidCurvatureModel, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontExactEllipsoidCurvatureModel, Dictionary);

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //
scalar frontExactEllipsoidCurvatureModel::devCosine(const scalar& longitude, const scalar& latitude, const point& p) const
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

bool frontExactEllipsoidCurvatureModel::converged(const scalar& longitude, const scalar& latitude, const point& p) const
{
    auto dev = devCosine(longitude, latitude, p);
    auto dist = distance(longitude, latitude, p);
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

point frontExactEllipsoidCurvatureModel::ellipsoidPoint(const scalar& longitude, const scalar& latitude) const
{
    return vector{
                    semiAxes_.x()*cos(longitude)*sin(latitude),
                    semiAxes_.y()*sin(longitude)*sin(latitude),
                    semiAxes_.z()*cos(latitude)
                 };
}

scalar frontExactEllipsoidCurvatureModel::distance(const scalar& longitude, const scalar& latitude, const point& p) const
{
    return mag(ellipsoidPoint(longitude, latitude) - p);
}

scalar frontExactEllipsoidCurvatureModel::curvature(const scalar& longitude, const scalar& latitude) const
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

frontExactEllipsoidCurvatureModel::interval frontExactEllipsoidCurvatureModel::findParameters(const point& p) const
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
                currentDistance = distance(u, v, p);

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontExactEllipsoidCurvatureModel::frontExactEllipsoidCurvatureModel(const dictionary& configDict)
:
    frontExactCurvatureModel(configDict),
    semiAxes_{configDict.lookup("semiAxes")},
    centre_{configDict.lookup("centre")},
    tolerance_{readScalar(configDict.lookup("tolerance"))}
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //
scalar frontExactEllipsoidCurvatureModel::curvatureAtPoint(const point& p) const
{
    // Move coordinate system to ellipse centre and project point to 
    // positive quadrant --> simplifies dealing with the trigonometric functions
    // This approach exploits symmetry
    point q = p - centre_;
    q.x() = mag(q.x());
    q.y() = mag(q.y());
    q.z() = mag(q.z());

    auto parameters = findParameters(q);
    auto longitude = parameters[0];
    auto latitude = parameters[1];

    return curvature(longitude, latitude);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
