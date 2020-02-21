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
    Foam::simpleTriangulator

SourceFiles
    simpleTriangulator.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Triangulates a given point set by decomposition using the geometric
    centre of the point set.
    
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

#ifndef simpleTriangulator_H
#define simpleTriangulator_H

#include <utility>

#include "fvCFD.H"

#include "analyticalPlane.H"
#include "analyticalSurface.H"

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
    explicit simpleTriangulator(pointField& vertices, List<triFace>& triangles);


    //- Destructor
    virtual ~simpleTriangulator();


    // Member Functions
    void setVertices(pointField& vertices);
    void setTriangleList(List<triFace>& triangles);
    void triangulateFace(labelList& facePointIDs, const vector& faceNormal);
    void triangulateFace(labelList& facePointIDs, const tmp<analyticalSurface>& surfaceTmp);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
