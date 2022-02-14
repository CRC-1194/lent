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
    Foam::frontConstructor

SourceFiles
    frontConstructor.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Constructs a triSurface by triangulating a given analytical surface on a
    per-cell basis with the given volume mesh.

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

#ifndef frontConstructor_H
#define frontConstructor_H

#include "fvCFD.H"
#include "pointMesh.H"
#include "triSurface.H"

#include "analyticalPlane.H"
#include "analyticalSurface.H"
#include "simpleTriangulator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontConstructor Declaration
\*---------------------------------------------------------------------------*/

class frontConstructor
{
    // Private data
    List<triFace> frontTriangles_;
    labelList intersectedEdges_;
    labelList intersectedCells_;
    labelList triaToCell_;
    pointField frontVertices_;

    const tmp<analyticalSurface> surfaceTmp_;
    const fvMesh& mesh_;


    // Private Member Functions
    bool differentSign(scalar, scalar) const;
    bool notInList(label, labelList&) const;
    label pointIndex(label, const labelList&) const;

    void findIntersectedEdges();
    void findIntersectedCells();
    void computeIntersections();
    void cellIntersections(const labelList&, labelList&) const;
    void createTriangles();


public:
    
    // Constructors
    explicit frontConstructor(const tmp<analyticalSurface>&, const fvMesh&);


    //- Destructor
    ~frontConstructor();


    // Member Functions
    const labelList& intersectedEdges() const
    {
        return intersectedEdges_;
    }

    const labelList& intersectedCells() const
    {
        return intersectedCells_;
    }

    const labelList& triangleToCell() const
    {
        return triaToCell_;
    }

    auto surface() const
    {
        return surfaceTmp_;
    }

    triSurface createTriSurface()
    {
        return triSurface(frontTriangles_, frontVertices_);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
