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
    Foam::FrontTracking:analyticalPlane

SourceFiles
    analyticalPlane.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for a plane.

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

#include "analyticalPlane.H"
#include "addToRunTimeSelectionTable.H"

#include <iomanip>

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalPlane, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalPlane, Dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
vector analyticalPlane::normalize(const vector& normalVector) const
{
    if (mag(normalVector) > SMALL)
    {
        return normalVector / mag(normalVector);
    }
    else
    {
        return normalVector / (mag(normalVector) + SMALL);
    }
}

void analyticalPlane::updateDistanceToOrigin()
{
    distanceOrigin_ = unitNormal_ & refPoint_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalPlane::analyticalPlane(const dictionary& configDict)
:
    analyticalSurface{configDict},
    refPoint_{configDict.get<point>("referencePoint")},
    unitNormal_{configDict.get<vector>("normalVector")}
{
    name_ = configDict.lookupOrDefault<word>("name", analyticalPlane::typeName);
    unitNormal_ = normalize(unitNormal_);
    updateDistanceToOrigin();
}

analyticalPlane::analyticalPlane(
    const point& refPoint,
    const vector& normal,
    const word name
)
:
    analyticalSurface{},
    refPoint_{refPoint},
    unitNormal_{normal},
    name_{name}
{
    unitNormal_ = normalize(unitNormal_);
    updateDistanceToOrigin();
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
    return (trialPoint - ((trialPoint - refPoint_) & (unitNormal_*unitNormal_)));
}

vector analyticalPlane::normalToPoint(const point&) const
{
    return unitNormal_;
}

point analyticalPlane::intersection(const point& pointA, const point& pointB) const
{
        point intersect(0.0, 0.0, 0.0);

        scalar distanceRatio = distance(pointA) / (distance(pointA)
                                 + distance(pointB));

        // In essence, use distance weighted central differencing scheme to
        // find intersetion with plane
        intersect = distanceRatio*pointB + (1.0 - distanceRatio) * pointA;

        return intersect;
}

void analyticalPlane::normal(const vector& newNormal)
{
    unitNormal_ = normalize(newNormal);
    updateDistanceToOrigin();
}

void analyticalPlane::referencePoint(const point& newRefPoint)
{
    refPoint_ = newRefPoint;
    updateDistanceToOrigin();
}

void analyticalPlane::writeParameters(const word fileName) const
{
    auto outputFile = outputStream(fileName);

    outputFile.stdStream() << std::setprecision(15);

    outputFile << "-------------------------------\n"
               << "type " << this->type() << '\n'
               << "referencePoint " << refPoint_ << '\n'
               << "normal " << unitNormal_ << '\n';
}


// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //
analyticalPlane& analyticalPlane::operator=(const analyticalPlane& rhs) 
{
    if (this != &rhs)
    {
        refPoint_ = rhs.refPoint_;
        unitNormal_ = rhs.unitNormal_;
    }

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
