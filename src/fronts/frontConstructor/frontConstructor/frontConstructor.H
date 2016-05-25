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
    Foam::frontConstructor

SourceFiles
    frontConstructor.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Constructs a triSurface by triangulating a given analytical surface on a
    per-cell basis with the given volume mesh.
    
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
    frontConstructor() = default;
    explicit frontConstructor(const tmp<analyticalSurface>, const fvMesh&);


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
