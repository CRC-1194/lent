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

#include "frontExactEllipseCurvatureModel.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontExactEllipseCurvatureModel, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontExactEllipseCurvatureModel, Dictionary);

    using constant::mathematical::pi;

// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //
point frontExactEllipseCurvatureModel::ellipsePoint(const scalar& angle) const
{
    return point{xSemiAxis_*cos(angle), ySemiAxis_*sin(angle), 0.0};
}

scalar frontExactEllipseCurvatureModel::distance(const scalar& angle, const point& p) const
{
    return mag(p - ellipsePoint(angle));
}

scalar frontExactEllipseCurvatureModel::curvature(const scalar& angle) const
{
    auto a = xSemiAxis_*sin(angle);
    auto b = ySemiAxis_*cos(angle);
    auto root = sqrt(a*a + b*b);

    return -xSemiAxis_*ySemiAxis_ / (root*root*root);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontExactEllipseCurvatureModel::frontExactEllipseCurvatureModel(const dictionary& configDict)
:
    frontExactCurvatureModel(configDict),
    xSemiAxis_{readScalar(configDict.lookup("xSemiAxis"))},
    ySemiAxis_{readScalar(configDict.lookup("ySemiAxis"))},
    centre_{configDict.lookup("centre")}
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //
scalar frontExactEllipseCurvatureModel::curvatureAtPoint(const point& p) const
{
    // Move coordinate system to ellipse centre and project point to 
    // positive quadrant --> simplifies dealing with the trigonometric functions
    // This approach exploits symmetry
    point q = p - centre_;
    q.x() = mag(q.x());
    q.y() = mag(q.y());
    q.z() = 0.0;

    scalar previousDistance = GREAT;
    scalar currentDistance = 0.0;
    scalar angle = 0.0;
    scalar increment = pi/20000.0;

    while (angle <= pi/2.0)
    {
        currentDistance = distance(angle, q);

        if (currentDistance > previousDistance)
        {
            angle -= increment;
            break;
        }

        previousDistance = currentDistance;
        angle += increment;
    }

    return curvature(angle);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
