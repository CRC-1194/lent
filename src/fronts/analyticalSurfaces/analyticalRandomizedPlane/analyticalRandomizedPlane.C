/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
    \\  /    A nd           | Copyright held by original author
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
    Foam::analyticalRandomizedPlane

SourceFiles
    analyticalRandomizedPlane.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    An analytical plane for which the normal and the reference point
    are computed from reference values plus random values

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

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
