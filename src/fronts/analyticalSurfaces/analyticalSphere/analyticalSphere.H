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
    Foam::FrontTracking:analyticalSphere

SourceFiles
    analyticalSphere.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for a sphere.

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

#ifndef analyticalSphere_H
#define analyticalSphere_H

#include <cmath>

#include "analyticalSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalSphere Declaration
\*---------------------------------------------------------------------------*/

class analyticalSphere
:
    public analyticalSurface
{
    // Private data
    point centre_;
    scalar radius_;
    word name_; 

public:

    TypeName ("sphere");

    // Constructors
    analyticalSphere() = default;
    analyticalSphere(const dictionary& configDict);
    analyticalSphere(
        const point& centre, 
        const scalar radius, 
        const word name = analyticalSurface::typeName
    );

    //- Destructor
    virtual ~analyticalSphere() {};


    // Member Functions
    virtual scalar distance(const point& trialpoint) const;
    virtual scalar signedDistance(const point& trialPoint) const;

    inline scalar value(const point& trialPoint) const
    {
        return signedDistance(trialPoint); 
    }

    inline vector gradient(const point& trialPoint) const
    {
        return normalToPoint(trialPoint); 
    }

    word name() const
    {
        return name_; 
    }

    virtual point normalProjectionToSurface(point& trialPoint) const;
    virtual vector normalToPoint(const point& trialPoint) const;
    virtual point intersection(const point& pointA, const point& pointB) const;
    virtual scalar curvatureAt(const point&) const
    {
        return -2.0/radius();
    }

    virtual point centre() const
    {
        return centre_;
    }

    virtual scalar radius() const
    {
        return radius_;
    }

    virtual void centre(const point& newCentre);
    virtual void radius(const scalar newRadius);

    virtual void writeParameters(const word fileName) const;

    // Member operators
    virtual analyticalSphere& operator=(const analyticalSphere& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
