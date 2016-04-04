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
    Foam::frontConstructor

Description

SourceFiles
    frontConstructor.C

\*---------------------------------------------------------------------------*/

#ifndef frontConstructor_H
#define frontConstructor_H

#include "fvCFD.H"

#include "analyticalSurface.H"
#include "analyticalPlane.H"
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
        pointField frontVertices_;

        const tmp<analyticalSurface> surfaceTmp_;
        const fvMesh& mesh_;


    // Private Member Functions
        bool differentSign(scalar, scalar) const;
        bool notInList(label, labelList&) const;
        label pointIndex(label, const labelList&) const;
        void signedDistance();

        void findIntersectedEdges();
        void findIntersectedCells();
        void computeIntersections();
        void cellIntersections(const labelList&, labelList&) const;
        void createTriangles();


public:
    
    // Constructors

        //- Construct null
        frontConstructor() = default;

        //- Construct from components
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

        // TODO: figure out a sensible and efficient way to contruct
        // a triSurfaceFront and pass it
        // FIXME: following functions are only temporary for testing purposes
        void workWork()
        {
            createTriangles();
        }

        const List<triFace>& triangles()
        {
            return frontTriangles_;
        }

        const pointField& vertices()
        {
            return frontVertices_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
