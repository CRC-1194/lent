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
    Foam::analyticalRandomizedEllipsoid

SourceFiles
    analyticalRandomizedEllipsoid.C

Author
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    Analytical Ellipsoiid with randomized parameters.

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
    centrePerturbation_{configDict.lookup("centrePerturbation")},
    semiAxesPerturbation_{configDict.lookup("semiAxesPerturbation")}
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

    forAll(randomSemiAxes, I)
    {
        if (randomSemiAxes[I] < SMALL)
        {
            randomSemiAxes[I] = mag(randomSemiAxes[I]) + SMALL;
        }
    }

    centre(randomCentre);
    semiAxes(randomSemiAxes);
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
