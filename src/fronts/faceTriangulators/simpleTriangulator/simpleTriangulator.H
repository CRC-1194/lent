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
    Foam::simpleTriangulator

Description

SourceFiles
    simpleTriangulator.C

\*---------------------------------------------------------------------------*/

#ifndef simpleTriangulator_H
#define simpleTriangulator_H

#include "fvCFD.H"

#include "analyticalSurface.H"
#include "analyticalPlane.H"

#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class simpleTriangulator Declaration
\*---------------------------------------------------------------------------*/

class simpleTriangulator
{
    // Private data
        pointField& vertices_;
        List<triFace>& triangles_;


    // Private Member Functions
        point geometricCentre(const labelList& pointIDs) const;
        scalar angle(const vector& refEdge, const vector& a,
                     scalar signedDistanceA) const;
        void linkedSort(scalarList& reference, labelList& dependent) const;
        void orderPoints(labelList& pointIDs, const point& refPoint,
                         const vector& normal) const;
        scalar unitLimiter(scalar) const;

        // Implements the actual triangulation algorithm. For now,
        // use a decomposition at the geometric centre
        virtual void triangulate(labelList& pointIDs, const vector& faceNormal);


public:

    // Constructors

        //- Construct null
        simpleTriangulator() = default;

        //- Construct from components
        explicit simpleTriangulator(pointField& vertices, List<triFace>& triangles);


    //- Destructor
    virtual ~simpleTriangulator();


    // Member Functions
        
        // Access
        void setVertices(pointField& vertices);
        void setTriangleList(List<triFace>& triangles);

        // Operation
        void triangulateFace(labelList& facePointIDs, const vector& faceNormal);
        void triangulateFace(labelList& facePointIDs, 
                             tmp<analyticalSurface> surfaceTmp);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
