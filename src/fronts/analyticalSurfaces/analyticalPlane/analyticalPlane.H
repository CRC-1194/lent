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

#ifndef analyticalPlane_H
#define analyticalPlane_H

#include "analyticalSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalPlane Declaration
\*---------------------------------------------------------------------------*/

class analyticalPlane
:
    public analyticalSurface
{
    // Private data
    point refPoint_;
    vector unitNormal_;
    scalar distanceOrigin_;
    word name_;

protected:
    vector normalize(const vector& normalVector) const;
    void updateDistanceToOrigin();

public:

    TypeName ("plane");

    // Constructors
    analyticalPlane() = default;
    analyticalPlane(const dictionary& configDict);
    analyticalPlane(
            const point& refPoint,
            const vector& normal,
            const word name = analyticalPlane::typeName
    );

    //- Destructor
    virtual ~analyticalPlane() {};


    // Member Functions
    virtual scalar distance(const point& trialPoint) const;
    virtual scalar signedDistance(const point& trialPoint) const;
    virtual point normalProjectionToSurface(point& trialPoint) const;
    virtual vector normalToPoint(const point& trialPoint) const;
    virtual point intersection(const point& pointA, const point& pointB) const;
    virtual scalar curvatureAt(const point&) const
    {
        return 0.0;
    }

    virtual vector normal() const
    {
        return unitNormal_;
    }

    virtual point referencePoint() const
    {
        return refPoint_;
    }

    inline scalar value(const point& trialPoint) const
    {
        return signedDistance(trialPoint); 
    }


    word name() const
    {
        return name_; 
    }

    virtual void normal(const vector& newNormal);
    virtual void referencePoint(const point& newRefPoint);

    virtual void writeParameters(const word fileName) const; 

    // Member operators
    virtual analyticalPlane& operator=(const analyticalPlane& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
