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
    Foam::analyticalSurface

SourceFiles
    analyticalSurface.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Abstract base class for definition of (simple) analytical surfaces.
    Intended for directly constructing a front or setting an exact signed
    distance field.

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

#ifndef analyticalSurface_H
#define analyticalSurface_H

#include "fvMesh.H"
#include "refCount.H"
#include "tmp.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking{

/*---------------------------------------------------------------------------*\
                         Class analyticalSurface Declaration
\*---------------------------------------------------------------------------*/

class analyticalSurface
:
    public refCount 
{

public:

    TypeName ("abstract");

    // Constructors
    analyticalSurface() = default;
    explicit analyticalSurface(const dictionary& configDict);

    declareRunTimeSelectionTable (
        tmp,
        analyticalSurface,
        Dictionary,
        (const dictionary& configDict),
        (configDict)
    );

    // Selectors
    static tmp<analyticalSurface> New(const dictionary& configDict);

    // Destructor
    virtual ~analyticalSurface() = default;

    // Member functions
    virtual scalar distance(const point&) const = 0;
    virtual scalar signedDistance(const point&) const = 0;
    virtual point normalProjectionToSurface(point&) const = 0;
    virtual vector normalToPoint(const point&) const = 0;
    virtual point intersection(const point&, const point&) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //