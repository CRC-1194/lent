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
    Foam::FrontTracking:analyticalRandomizedSphere

SourceFiles
    analyticalRandomizedSphere.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for a randomized sphere.

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

#include "analyticalRandomizedSphere.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalRandomizedSphere, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalRandomizedSphere, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalRandomizedSphere::analyticalRandomizedSphere(const dictionary& configDict)
:
    analyticalSphere{configDict},
    centrePerturbation_{configDict.get<vector>("centrePerturbation")},
    radiusPerturbation_{configDict.get<scalar>("radiusPerturbation")}
{
    originalCentre_ = centre();
    originalRadius_ = radius();

    randomize();
}

analyticalRandomizedSphere::analyticalRandomizedSphere(const point& centre, const scalar radius, const point& centrePerturbation, const scalar radiusPerturbation)
:
   analyticalSphere{centre, radius},
   originalCentre_{centre},
   originalRadius_{radius},
   centrePerturbation_{centrePerturbation},
   radiusPerturbation_{radiusPerturbation}
{
    randomize();
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void analyticalRandomizedSphere::randomize()
{
    auto randomCentre = originalCentre_ + noiseGen_.noise<vector>(centrePerturbation_);
    auto randomRadius = originalRadius_ + noiseGen_.noise<scalar>(radiusPerturbation_);

    centre(randomCentre);
    
    if (randomRadius > SMALL)
    {
        radius(randomRadius);
    }
}

void analyticalRandomizedSphere::randomPlacementIn(const boundBox& bbox)
{
    auto maxDisplacement = 0.5*bbox.span() - originalRadius_*vector{1,1,1};
    auto newCentre = bbox.midpoint() + noiseGen_.noise<vector>(maxDisplacement);

    centre(newCentre);
}


// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //
analyticalRandomizedSphere& analyticalRandomizedSphere::operator=(const analyticalRandomizedSphere& rhs)
{
    if (this != &rhs)
    {
        analyticalSphere::operator=(rhs);
        originalCentre_ = rhs.originalCentre_;
        originalRadius_ = rhs.originalRadius_;
        centrePerturbation_ = rhs.centrePerturbation_;
        radiusPerturbation_ = rhs.radiusPerturbation_;
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
