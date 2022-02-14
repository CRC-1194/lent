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
    Foam::edgeInversedNormalCalculator

SourceFiles
    edgeInversedNormalCalculator.C

Authors:
    Tobias Tolle   (tolle@mma.tu-darmstadt.de)

Description

    Compute the normals at the front vertices in two steps. First, the normals
    on the edges are approximated using an average of the two connected
    triangles. Each normal is assigned the area of the opposite triangle as
    weight. Second, the vertex normal is computed by averaging the edge
    normals of the connected edges. The inversed edge lengths (inversed
    distance weighting) are used as weights.

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

#ifndef edgeInversedNormalCalculator_H
#define edgeInversedNormalCalculator_H

#include "frontVertexNormalCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class edgeInversedNormalCalculator Declaration
\*---------------------------------------------------------------------------*/

class edgeInversedNormalCalculator
:
    public frontVertexNormalCalculator
{

public:

    TypeName ("edgeNormalInversed");

    // Constructors
    edgeInversedNormalCalculator() = default;
    explicit edgeInversedNormalCalculator(const dictionary& configDict);

    //- Destructor
    virtual ~edgeInversedNormalCalculator() = default;

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
