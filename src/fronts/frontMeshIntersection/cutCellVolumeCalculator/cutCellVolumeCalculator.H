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
    Foam::cutCellVolumeCalculator

Description

SourceFiles
    cutCellVolumeCalculator.C

\*---------------------------------------------------------------------------*/

#ifndef cutCellVolumeCalculator_H
#define cutCellVolumeCalculator_H

#include "pointFields.H"

#include "lentCommunication.H"
#include "faceTriangulator.H"

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

public:

    // Constructors
    cutCellVolumeCalculator() = default;
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
