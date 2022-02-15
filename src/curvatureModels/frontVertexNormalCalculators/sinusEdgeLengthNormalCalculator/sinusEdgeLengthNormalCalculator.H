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
    Foam::sinusEdgeLengthNormalCalculator

SourceFiles
    sinusEdgeLengthNormalCalculator.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)

Description

    Compute the vertex normals as a weighted average of the normals of the
    connected triangles.
    The weights are computed as in interface tracking using as weight
    w = sin(alpha) / (|d0| * |d1|), where alpha is the angle at the vertex and
    d0 and d1 are the edges connecting the vertex with the two other points
    of the given triangle.

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
