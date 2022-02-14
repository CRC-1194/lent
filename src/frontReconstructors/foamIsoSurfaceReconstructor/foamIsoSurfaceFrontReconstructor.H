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
    Foam::foamIsoSurfaceFrontReconstructor

SourceFiles
    foamIsoSurfaceFrontReconstructor.C

Author
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Reconstruct a front from as an iso-surface. 

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


#ifndef foamIsoSurfaceFrontReconstructor_H
#define foamIsoSurfaceFrontReconstructor_H

#include "frontReconstructor.H"
#include "normalConsistency.H"
#include "frontSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class foamIsoSurfaceFrontReconstructor Declaration
\*---------------------------------------------------------------------------*/

class foamIsoSurfaceFrontReconstructor
    :
        public frontReconstructor
{

    scalar mergeTolerance_;
    Switch regularize_;

    tmp<normalConsistency> consistencyAlgTmp_;

public:

    TypeName ("foamIsoSurface");

    // Constructors
    foamIsoSurfaceFrontReconstructor(const dictionary& configDict);

    // Member Functions
    virtual void reconstructFront(
        triSurfaceFront& front,
        const volScalarField& signedDistance,
        const pointScalarField& pointSignedDistance
    ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
