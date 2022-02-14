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
    Foam::FrontTracking:analyticalRandomizedEllipsoid

SourceFiles
    analyticalRandomizedEllipsoid.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Specialization of the analyticalSurface class for a randomized ellipsoid.

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

#include "analyticalRandomizedEllipsoid.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalRandomizedEllipsoid, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalRandomizedEllipsoid, Dictionary);
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalRandomizedEllipsoid::analyticalRandomizedEllipsoid(const dictionary& configDict)
:
    analyticalEllipsoid{configDict},
    centrePerturbation_{configDict.get<vector>("centrePerturbation")},
    semiAxesPerturbation_{configDict.get<vector>("semiAxesPerturbation")}
{
    originalCentre_ = centre();
    originalSemiAxes_ = semiAxes();

    randomize();
}

analyticalRandomizedEllipsoid::analyticalRandomizedEllipsoid(const point& centre, const vector& semiAxes, const point& centrePerturbation, const vector& semiAxesPerturbation)
:
   analyticalEllipsoid{centre, semiAxes},
   originalCentre_{centre},
   originalSemiAxes_{semiAxes},
   centrePerturbation_{centrePerturbation},
   semiAxesPerturbation_{semiAxesPerturbation}
{
    randomize();
} 


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void analyticalRandomizedEllipsoid::randomize()
{
    auto randomCentre = originalCentre_ + noiseGen_.noise<vector>(centrePerturbation_);
    auto randomSemiAxes = originalSemiAxes_ + noiseGen_.noise<vector>(semiAxesPerturbation_);

    for (direction I = 0; I != 3; ++I)
    {
        if (randomSemiAxes[I] < SMALL)
        {
            randomSemiAxes[I] = mag(randomSemiAxes[I]) + SMALL;
        }
    }

    centre(randomCentre);
    semiAxes(randomSemiAxes);
}

void analyticalRandomizedEllipsoid::randomPlacementIn(const boundBox& bbox)
{
    auto maxDisplacement = 0.5*bbox.span() - originalSemiAxes_;
    auto newCentre = bbox.midpoint() + noiseGen_.noise<vector>(maxDisplacement);

    centre(newCentre);
}


// * * * * * * * * * * * * * * Member Operators    * * * * * * * * * * * * * * //
analyticalRandomizedEllipsoid& analyticalRandomizedEllipsoid::operator=(const analyticalRandomizedEllipsoid& rhs)
{
    if (this != &rhs)
    {
        analyticalEllipsoid::operator=(rhs);
        originalCentre_ = rhs.originalCentre_;
        originalSemiAxes_ = rhs.originalSemiAxes_;
        centrePerturbation_ = rhs.centrePerturbation_;
        semiAxesPerturbation_ = rhs.semiAxesPerturbation_;
    }

    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
