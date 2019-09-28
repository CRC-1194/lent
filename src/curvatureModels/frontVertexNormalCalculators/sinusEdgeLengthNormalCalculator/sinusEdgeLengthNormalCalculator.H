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
    Foam::sinusEdgeLengthNormalCalculator

SourceFiles
    sinusEdgeLengthNormalCalculator.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Compute the vertex normals as a weighted average of the normals of the
    connected triangles.
    The weights are computed as in interface tracking using as weight
    w = sin(alpha) / (|d0| * |d1|), where alpha is the angle at the vertex and
    d0 and d1 are the edges connecting the vertex with the two other points
    of the given triangle.

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

#ifndef sinusEdgeLengthNormalCalculator_H
#define sinusEdgeLengthNormalCalculator_H

#include "frontVertexNormalCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class sinusEdgeLengthNormalCalculator Declaration
\*---------------------------------------------------------------------------*/

class sinusEdgeLengthNormalCalculator
:
    public frontVertexNormalCalculator
{

private:
    FixedList<label,2> notCentrePoint(const label& centreID, const face& trianglePointLabels) const;

public:

    TypeName ("sinusEdgeLength");

    // Constructors
    sinusEdgeLengthNormalCalculator() = default;
    explicit sinusEdgeLengthNormalCalculator(const dictionary& configDict);

    //- Destructor
    virtual ~sinusEdgeLengthNormalCalculator() = default;

    // Member Functions
    virtual tmp<triSurfaceFrontPointVectorField> vertexNormals(const fvMesh&, const triSurfaceFront&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
