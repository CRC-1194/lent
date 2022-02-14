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
    Foam::maxNormalAngleFrontReconstructionModel

SourceFiles
    maxNormalAngleFrontReconstructionModel.C

Authors:
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Invoke front reconstruction if the maximum area of a front triangle
    exceeds a given threshold.
    Based on the cell volume of the Eulerian Mesh a characteristic area
    is computed and scaled with user prescribed factor.

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

#ifndef maxAreaFrontReconstructionModel_H
#define maxAreaFrontReconstructionModel_H

#include "frontReconstructionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: extend this class so that also a minimum area can be prescribed
// and checked. First check if this makes sense
namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class maxAreaFrontReconstructionModel Declaration
\*---------------------------------------------------------------------------*/

class maxAreaFrontReconstructionModel
:
    public frontReconstructionModel
{
  
    mutable scalar maxArea_;
    scalar scaleFactor_;

    void computeMaxArea(const fvMesh& mesh) const;

public:

    TypeName ("maxArea");

    // Constructors
    explicit maxAreaFrontReconstructionModel(const dictionary& configDict);
    

    // Member Functions
    virtual bool reconstructionRequired(
        const triSurfaceFront& front,
        const volScalarField& signedDistance
    ) const;

    scalar maxArea() const
    {
        return maxArea_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
