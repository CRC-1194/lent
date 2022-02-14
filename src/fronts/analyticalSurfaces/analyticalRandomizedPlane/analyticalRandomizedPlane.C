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

#include "analyticalRandomizedPlane.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(analyticalRandomizedPlane, 0);
    addToRunTimeSelectionTable(analyticalSurface, analyticalRandomizedPlane, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
analyticalRandomizedPlane::analyticalRandomizedPlane(const dictionary& configDict)
:
    analyticalPlane{configDict},
    refPointPerturbation_{configDict.get<point>("referencePointPerturbation")},
    normalPerturbation_{configDict.get<vector>("normalVectorPerturbation")}
{
    originalRefPoint_ = referencePoint();
    originalNormal_ = normal();

    randomize();
}

analyticalRandomizedPlane::analyticalRandomizedPlane(const point& refPoint, const vector& normal, const point& refPointPerturbation, const vector& normalPerturbation)
:
    analyticalPlane{refPoint, normal},
    originalRefPoint_{refPoint},
    originalNormal_{refPoint},
    refPointPerturbation_{refPointPerturbation},
    normalPerturbation_{normalPerturbation}
{
    originalNormal_ = normalize(originalNormal_);
    randomize();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void analyticalRandomizedPlane::randomize()
{
    auto randomRefPoint = originalRefPoint_ + noiseGen_.noise<vector>(refPointPerturbation_);
    auto randomNormal = originalNormal_ + noiseGen_.noise<vector>(normalPerturbation_);
    randomNormal = normalize(randomNormal);

    referencePoint(randomRefPoint);
    normal(randomNormal);
}

void analyticalRandomizedPlane::randomPlacementIn(const boundBox& bbox)
{
    referencePoint(bbox.midpoint() + noiseGen_.noise<vector>(0.5*bbox.span()));
}

// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //
analyticalRandomizedPlane& analyticalRandomizedPlane::operator=(const analyticalRandomizedPlane& rhs) 
{
    if (this != &rhs)
    {
        originalRefPoint_ = rhs.originalRefPoint_;
        originalNormal_ = rhs.originalNormal_;
        refPointPerturbation_ = rhs.refPointPerturbation_;
        normalPerturbation_ = rhs.normalPerturbation_;
    }

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
