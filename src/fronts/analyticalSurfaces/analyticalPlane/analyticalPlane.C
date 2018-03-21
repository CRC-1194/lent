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
    Foam::analyticalPlane

SourceFiles
    analyticalPlane.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Specialization of the analyticalSurface class for a plane.

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalPlane::analyticalPlane(const dictionary& configDict)
:
    analyticalSurface{configDict},
    refPoint_{configDict.lookup("referencePoint")},
    unitNormal_{configDict.lookup("normalVector")}
{
    unitNormal_ = normalize(unitNormal_);
    distanceOrigin_ = unitNormal_ & refPoint_;
}

analyticalPlane::analyticalPlane(const point& refPoint, const vector& normal)
:
    analyticalSurface{},
    refPoint_{refPoint},
    unitNormal_{normal}
{
    unitNormal_ = normalize(unitNormal_);
    distanceOrigin_ = unitNormal_ & refPoint_;
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
    point projected(0.0, 0.0, 0.0);

    // Projection method fails for zero vector, thus provide 
    // alternative method for this case
    if (fabs(trialPoint & trialPoint) > SMALL)
    {
        scalar projectedDistanceToOrigin = unitNormal_ & trialPoint;

        projected = projectedDistanceToOrigin / distanceOrigin_ * trialPoint;
    }
    else
    {
        projected = distanceOrigin_ * unitNormal_;
    }

    return projected;
}

vector analyticalPlane::normalToPoint(const point& trialPoint) const
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
    unitNormal_ = newNormal / (mag(newNormal) + SMALL);
}

void analyticalPlane::referencePoint(const point& newRefPoint)
{
    refPoint_ = newRefPoint;
}

void analyticalPlane::writeParameters(const word fileName) const
{
    // TODO: when porting to the latest OpenFOAM-plus version,
    // open file in append mode so the entire history of surface
    // parameters is saved (TT)
    OFstream outputFile(fileName);

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
