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
    Foam::convexFrontNormalConsistency

SourceFiles
    convexFrontNormalConsistency.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Ensures normal consistency by computing a geometric centre.
    ONLY applicable to CONVEX fronts.

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


#ifndef convexFrontNormalConsistency_H
#define convexFrontNormalConsistency_H

#include "normalConsistency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class convexFrontNormalConsistency Declaration
\*---------------------------------------------------------------------------*/

class convexFrontNormalConsistency
    :
        public normalConsistency
{

private:
    
    word orientation_;
    scalar orientationSign_;
    
    // Private Member Functions
    virtual void runNormalConsistencyAlgorithm(
        triSurfaceFront& front,
        const volScalarField& signedDistance,       // Not used by this alg.
        const pointScalarField& pointSignedDistance // Not used by this alg.
    ) const;

public:

    TypeName ("convexFrontNormal");

    // Constructors
        convexFrontNormalConsistency(const dictionary& configDict);

    // Selectors

        static tmp<convexFrontNormalConsistency> New(const dictionary& configDict);

    // Destructor
        virtual ~convexFrontNormalConsistency() {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
