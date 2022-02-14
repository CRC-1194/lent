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
    Foam::analyticalSurfaceNormalConsistency

SourceFiles
    analyticalSurfaceNormalConsistency.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Uses an analytical surface description to ensure consistent normal
    orientation.
    Intended for testing purposes only.

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


#ifndef analyticalSurfaceNormalConsistency_H
#define analyticalSurfaceNormalConsistency_H

#include "analyticalSurface.H"
#include "normalConsistency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalSurfaceNormalConsistency Declaration
\*---------------------------------------------------------------------------*/

class analyticalSurfaceNormalConsistency
    :
        public normalConsistency
{

private:
    
    bool reverseNormal_{false};
    tmp<analyticalSurface> surfaceTmp_;
    
    // Private Member Functions
    virtual void runNormalConsistencyAlgorithm(
        triSurfaceFront& front,
        const volScalarField& signedDistance,       // Not used by this alg.
        const pointScalarField& pointSignedDistance // Not used by this alg.
    ) const;

public:

    TypeName ("analyticalSurfaceNormal");

    // Constructors
        analyticalSurfaceNormalConsistency(const dictionary& configDict);
        analyticalSurfaceNormalConsistency(const analyticalSurface& surface, const bool reverseNormal=false);


    // Destructor
        virtual ~analyticalSurfaceNormalConsistency() {};
        
    // Member Functions
        bool makeFrontNormalsConsistent(triSurfaceFront& front) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
