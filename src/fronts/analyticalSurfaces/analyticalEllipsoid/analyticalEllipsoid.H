/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 Tobias Tolle, TU Darmstadt 
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

\*---------------------------------------------------------------------------*/

#ifndef analyticalEllipsoid_H
#define analyticalEllipsoid_H

#include <cmath>
#include <functional>

#include "analyticalSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalEllipsoid Declaration
\*---------------------------------------------------------------------------*/

class analyticalEllipsoid
:
    public analyticalSurface
{

public:
    using parameterPair = FixedList<scalar,2>;

private:
    
    // Private data
    point centre_;
    vector oneBySemiAxisSqr_;

    // NOTE: all private member functions functions assume that the given point
    // arguments are in the elliposid's reference frame. The axes orientation is the
    // same as in the global coordinate system, but the origin is located at the
    // centroid of the ellipsoid.
    scalar levelSetValueOf(const point& aPoint) const;
    vector levelSetGradientAt(const point& aPoint) const;
    point moveToReferenceFrame(const point& aPoint) const;
    parameterPair intersectEllipsoidWithLine(const point& refPoint, const vector& path) const;
    scalar ellipsoidCurvature(const point& p) const;
    label minorSemiAxisIndex() const;
    // TODO: extract bisection method since it may be used else wehere (TT)
    scalar bisection(const std::function<scalar(scalar)>& rootFunction, parameterPair interval) const;

    // Counterparts of public interface functions formulated in the
    // ellipsoid reference system
    scalar signedDistanceRefSytem(const point P) const;
    point normalProjectionToSurfaceRefSystem(point P) const;
    vector normalToPointRefSystem(const point P) const;
    point intersectionRefSystem(const point A, const point B) const;
    scalar curvatureAtRefSystem(const point P) const;
     

public:

    TypeName ("ellipsoid");

    // Constructors
    analyticalEllipsoid() = default;
    analyticalEllipsoid(const dictionary& configDict);
    analyticalEllipsoid(const point& centre, const vector& semiAxes);

    //- Destructor
    virtual ~analyticalEllipsoid() {};


    // Member Functions
    virtual scalar distance(const point& trialPoint) const;
    virtual scalar signedDistance(const point& trialPoint) const;
    virtual point normalProjectionToSurface(point& trialPoint) const;
    virtual vector normalToPoint(const point& trialPoint) const;
    virtual point intersection(const point& pointA, const point& pointB) const;
    virtual scalar curvatureAt(const point& p) const;

    virtual point centre() const;
    virtual vector semiAxes() const;
    
    // Set new parameters for the ellipsoid
    virtual void oneBySemiAxisSquare(const vector& oneBySemiAxisSquare);
    virtual void centre(const vector& newCentre);
    virtual void semiAxes(const vector& newSemiAxes);

    virtual void writeParameters(const word fileName) const;

    // Member operators
    virtual analyticalEllipsoid& operator=(const analyticalEllipsoid& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
