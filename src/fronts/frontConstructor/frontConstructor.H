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
    frontConstructorI.H
    frontConstructor.C
    frontConstructorIO.C

\*---------------------------------------------------------------------------*/

#ifndef frontConstructor_H
#define frontConstructor_H

#include "fvCFD.H"

#include "analyticalSurface.H"
#include "analyticalPlane.H"

#include <utility>

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
    // TODO: there has to be a mapping relating each cell to a list of 
    // contained triangles. Preferable way is to use the lentCommunication
    // somehow
    // TODO: add the faceTriangulation class as private data when available

    const tmp<analyticalSurface> surfaceTmp_;
    const fvMesh& mesh_;


    // Private Member Functions
    bool differentSign(scalar, scalar) const;
    bool notInList(label, labelList&) const;
    label pointIndex(label, labelList&) const;
    void createTriangles(const labelList&, const label&);
    void signedDistance();

    void findIntersectedEdges();
    void findIntersectedCells();
    void computeIntersections();

    // These belong to the to-be-created faceTriangulation class
    /*
    point geometricCentre(const labelList&, const pointField&) const;
    scalar angle(const vector&, const vector&, scalar) const;
    void linkedSort(scalarList&, labelList&) const;
    void orderPoints(labelList&, const pointField&, const point&,
                        const vector&) const;
    */

public:
    
    // Constructors

        //- Construct null
        frontConstructor() = default;

        //- Construct from components
        explicit frontConstructor(const tmp<analyticalSurface>, const fvMesh&);


    //- Destructor
    ~frontConstructor();


    // Member Functions
        const labelList& intersectedEdges()
        {
            return intersectedEdges_;
        }

        const labelList& intersectedCells()
        {
            return intersectedCells_;
        }



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
