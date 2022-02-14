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
    Foam::FrontTracking:analyticalRandomizedPlane

SourceFiles
    analyticalRandomizedPlane.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for a randomized plane.

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

#ifndef analyticalRandomizedPlane_H
#define analyticalRandomizedPlane_H

#include "noiseGenerator.H"

#include "analyticalPlane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalRandomizedPlane Declaration
\*---------------------------------------------------------------------------*/

class analyticalRandomizedPlane
:
    public analyticalPlane
{
    // Private data
    noiseGenerator noiseGen_;

    point originalRefPoint_;
    vector originalNormal_;
    point refPointPerturbation_;
    vector normalPerturbation_;

public:

    TypeName ("randomizedPlane");

    // Constructors
    analyticalRandomizedPlane() = default;
    analyticalRandomizedPlane(const dictionary& configDict);
    analyticalRandomizedPlane(const point& refPoint, const vector& normal, const point& refPointPerturbation, const vector& normalPerturbation);

    //- Destructor
    virtual ~analyticalRandomizedPlane() {};

    // Member Functions
    virtual void randomize();
    virtual void randomPlacementIn(const boundBox& bbox);

    // Member operators
    virtual analyticalRandomizedPlane& operator=(const analyticalRandomizedPlane& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
