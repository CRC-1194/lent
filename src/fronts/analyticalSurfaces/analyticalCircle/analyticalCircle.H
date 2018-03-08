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
    Foam::analyticalCircle

SourceFiles
    analyticalCircle.C

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

#ifndef analyticalCircle_H
#define analyticalCircle_H

#include <cmath>

#include "analyticalSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalCircle Declaration
\*---------------------------------------------------------------------------*/

class analyticalCircle
:
    public analyticalSurface
{
    // Private data
    point centre_;
    scalar radius_;
    tensor projector_;

public:

    TypeName ("circle");

    // Constructors
    analyticalCircle() = default;
    analyticalCircle(const dictionary& configDict);
    analyticalCircle(const point& centre, const scalar radius, const vector& emptyDirection);

    //- Destructor
    virtual ~analyticalCircle() {};


    // Member Functions
    virtual scalar distance(const point& trialpoint) const;
    virtual scalar signedDistance(const point& trialPoint) const;
    virtual point normalProjectionToSurface(point& trialPoint) const;
    virtual vector normalToPoint(const point& trialPoint) const;
    virtual point intersection(const point& pointA, const point& pointB) const;
    virtual scalar curvatureAt(const point& p) const
    {
        return -1.0/radius();
    }

    virtual point centre() const
    {
        return centre_;
    }

    virtual scalar radius() const
    {
        return radius_;
    }

    virtual void centre(const vector newCentre)
    {
        centre_ = newCentre;
    }

    virtual void radius(const scalar newRadius)
    {
        radius_ = newRadius;
    }

    // Member operators
    virtual analyticalCircle& operator=(const analyticalCircle& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
