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
    Foam::cutCellVolumeCalculator

SourceFiles
    cutCellVolumeCalculator.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Compute the phase volume in cells which are intersected by the front

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

#ifndef cutCellVolumeCalculator_H
#define cutCellVolumeCalculator_H

#include "pointFields.H"

#include "lentCommunication.H"
#include "simpleTriangulator.H"

#include <assert.h>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class cutCellVolumeCalculator Declaration
\*---------------------------------------------------------------------------*/

class cutCellVolumeCalculator
{
    // Private data
    const fvMesh& mesh_;
    const triSurfaceFront& front_;
    const word cellDistFieldName_;
    const word pointDistFieldName_;
    Map<List<label>> cellToTria_;
    Map<List<label>> cellToFace_;


    // Private Member Functions
    void addToMap(const label& cellID, const label& mappedID,
                  Map<List<label>>& map);
    void cellToTriangle();
    void cellToFace();
    scalar cutCellVolume(const label& cellID) const;
    scalar tetVolume(const pointField& points, const triFace& base,
                     const point& top) const;
    scalar polyhedraVolume(const pointField& points,
                           const List<triFace>& triangles) const;

    point geometricCentre(const pointField& points) const;
    void frontFragment(const labelList& frontTriangleIDs,
                       pointField& fragmentVertices,
                       List<triFace>& fragmentTriangles) const;
    label facePosition(const face& cellFace,
                       const pointScalarField& distance) const;
    label intersectionID(const point& a, const point& b,
                         const pointField& vertices, const label limit) const;
    bool vectorsParallel(const vector& a, const vector& b) const;
    vector provideNormal(const face& cellFace, const pointField& vertices) const;


public:

    // Constructors
    cutCellVolumeCalculator(const fvMesh& mesh, const triSurfaceFront& front,
                            const word& cellDistFieldName,
                            const word& pointDistFieldName);


    //- Destructor
    ~cutCellVolumeCalculator();


    // Member Functions
    scalar cellVolumePositivePhase(const label& cellIndex) const;
    scalar cellVolumeNegativePhase(const label& cellIndex) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
