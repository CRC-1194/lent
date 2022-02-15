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
    Foam::FrontTracking:analyticalRandomizedCircle

SourceFiles
    analyticalRandomizedCircle.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for a randomized circle.

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

#include "analyticalRandomizedCircle.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalRandomizedCircle, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalRandomizedCircle, Dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalRandomizedCircle::analyticalRandomizedCircle(const dictionary& configDict)
:
    analyticalCircle{configDict},
    centrePerturbation_{configDict.get<vector>("centrePerturbation")},
    radiusPerturbation_{configDict.get<scalar>("radiusPerturbation")}
{
    originalCentre_ = centre();
    originalRadius_ = radius();

    randomize();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void analyticalRandomizedCircle::randomize()
{
    auto perturbedCentre = originalCentre_ + noiseGen_.noise<vector>(centrePerturbation_);
    auto perturbedRadius = originalRadius_ + noiseGen_.noise<scalar>(radiusPerturbation_);
    centre(perturbedCentre);

    if (perturbedRadius > SMALL)
    {
        radius(perturbedRadius);
    }
}

void analyticalRandomizedCircle::randomPlacementIn(const boundBox& bbox)
{
    auto newCentre = bbox.midpoint() + 0.5*noiseGen_.noise<vector>(bbox.span() - originalRadius_*vector{1,1,1});

    centre(newCentre);
}


// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //
analyticalRandomizedCircle& analyticalRandomizedCircle::operator=(const analyticalRandomizedCircle& rhs)
{
    if (this != &rhs)
    {
        analyticalCircle::operator=(rhs);
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
