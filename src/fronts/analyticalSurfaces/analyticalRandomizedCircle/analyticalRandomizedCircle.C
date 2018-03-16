/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
    centrePerturbation_{configDict.lookup("centrePerturbation")},
    radiusPerturbation_{readScalar(configDict.lookup("radiusPerturbation"))}
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
