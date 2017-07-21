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
    Foam::maxNormalAngleFrontReconstructionModel

SourceFiles
    maxNormalAngleFrontReconstructionModel.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description
    Combines two front reconstruction models. Invokes reconstruction if
    either one of the criteria is fulfilled.

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

#include "addToRunTimeSelectionTable.H"

#include "twoCriteriaFrontReconstructionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(twoCriteriaFrontReconstructionModel, 0);
    addToRunTimeSelectionTable(frontReconstructionModel, twoCriteriaFrontReconstructionModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoCriteriaFrontReconstructionModel::twoCriteriaFrontReconstructionModel(const dictionary& configDict)
:
    frontReconstructionModel{configDict},
    firstCriterionTmp_{
        frontReconstructionModel::New(
            configDict.subDict("firstCriterion")
        )
    },
    secondCriterionTmp_{
        frontReconstructionModel::New(
            configDict.subDict("secondCriterion")
        )
    }
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool twoCriteriaFrontReconstructionModel::reconstructionRequired(
    const triSurfaceFront& front,
    const volScalarField& signedDistance
) const
{
    bool firstConditionFulfilled =
        firstCriterionTmp_->reconstructionRequired(front, signedDistance);

    bool secondConditionFulfilled =
        secondCriterionTmp_->reconstructionRequired(front, signedDistance);

    if (firstConditionFulfilled || secondConditionFulfilled)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
