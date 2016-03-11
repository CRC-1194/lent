/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::analyticalPlane

Description
    Specializes the analyticalSurface class for a plane

SourceFiles
    analyticalPlane.C

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

public:

    TypeName ("Plane");

    // Constructors
    analyticalPlane() = default;
    analyticalPlane(const dictionary& configDict);

    //- Destructor
    virtual ~analyticalPlane() {};


    // Member Functions
    virtual scalar distance(const point&) const;
    virtual scalar signedDistance(const point&) const;
    virtual point normalProjectionToSurface(point&) const;
    virtual vector normalToPoint(const point&) const;
    virtual point intersection(const point&, const point&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
