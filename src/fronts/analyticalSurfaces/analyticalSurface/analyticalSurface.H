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
    Foam::FrontTracking:analyticalSurface

SourceFiles
    analyticalSurfaceI.H
    analyticalSurface.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Abstract base class for definition of (simple) analytical surfaces.
    Intended for directly constructing a front or setting an exact signed
    distance field.

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

#ifndef analyticalSurface_H
#define analyticalSurface_H

#include "fvMesh.H"
#include "OFstream.H"
#include "refCount.H"
#include "tmp.H"
#include "triSurface.H"
#include "typeInfo.H"
#include "surfaceFields.H"
#include "volFields.H"


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
protected:
    OFstream outputStream(const word& fileName) const;

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
    )

    // Selectors
    static tmp<analyticalSurface> New(const dictionary& configDict);

    // Destructor
    virtual ~analyticalSurface() = default;

    // Member functions
    virtual scalar distance(const point&) const {return 0.0;};
    virtual scalar signedDistance(const point&) const {return 0.0;};

    template<class Container, class PointField>
    void setDistance(Container& distanceField, const PointField& points) const;

    void setDistance(volScalarField& signedDistance) const;
    void setDistance(surfaceScalarField& signedDistance) const;
    void setDistance(pointScalarField& signedDistance) const;

    virtual point normalProjectionToSurface(point&) const {return point{0,0,0};};
    virtual vector normalToPoint(const point&) const {return vector{0,0,0};};
    virtual point intersection(const point&, const point&) const {return point{0,0,0};};
    virtual void moveFrontToSurface(triSurface& front) const;
    virtual void makeNormalOrientationConsistent(triSurface& front, const bool outwardOrientation=true) const;
    virtual scalar curvatureAt(const point&) const {return 0.0;};

    virtual void randomize(){};
    virtual void randomPlacementIn(const boundBox&){};
    virtual void writeParameters(const word) const {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "analyticalSurfaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
