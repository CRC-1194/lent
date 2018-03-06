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
    Foam::analyticalEllipse

SourceFiles
    analyticalEllipse.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Specialization of the analyticalSurface class for an ellipse.

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

#ifndef analyticalEllipse_H
#define analyticalEllipse_H

#include "analyticalSurface.H"

#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalEllipse Declaration
\*---------------------------------------------------------------------------*/

class analyticalEllipse
:
    public analyticalSurface
{

public:
    using parameterPair = FixedList<scalar,2>;

private:

    // Private data
    point centre_;
    vector semiAxes_;
    tensor projector_;

    // Private member functions
    void ensureValidHalfAxes();
    void ensureValidCentre();

    label majorSemiAxisIndex() const;
    label minorSemiAxisIndex() const;

    scalar levelSetValueOf(const point& aPoint) const;
    vector levelSetGradientAt(const point& aPoint) const;

    point moveToReferenceFrame(const point& aPoint) const;
    vector emptyComponent(const point& aPoint) const;

    parameterPair intersectEllipseWithLine(const point& refPoint, const vector& path) const;
    scalar ellipseCurvature(const point& aPoint) const;

    // TODO: extract bisection method since it may be used else wehere (TT)
    scalar bisection(const std::function<scalar(scalar)>& rootFunction, parameterPair interval) const;

public:

    TypeName ("ellipse");

    // Constructors
    analyticalEllipse() = default;
    analyticalEllipse(const dictionary& configDict);
    analyticalEllipse(const point& centre, const vector& semiAxes, const vector& emptyDirection);

    //- Destructor
    virtual ~analyticalEllipse() {};


    // Member Functions
    virtual scalar distance(const point& trialPoint) const;
    virtual scalar signedDistance(const point& trialPoint) const;
    virtual point normalProjectionToSurface(point& trialPoint) const;
    virtual vector normalToPoint(const point& trialPoint) const;
    virtual point intersection(const point& pointA, const point& pointB) const;
    virtual scalar curvatureAt(const point& aPoint) const;

    virtual point centre() const
    {
        return centre_;
    }

    virtual const vector& semiAxes() const
    {
        return semiAxes_;
    }

    // Member operators
    virtual analyticalEllipse& operator=(const analyticalEllipse& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
