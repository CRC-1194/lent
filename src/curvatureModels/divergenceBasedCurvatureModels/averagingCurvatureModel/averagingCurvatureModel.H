/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 AUTHOR,AFFILIATION
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

Class
    Foam::averagingCurvatureModel

Description

SourceFiles
    averagingCurvatureModelI.H
    averagingCurvatureModel.C
    averagingCurvatureModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef averagingCurvatureModel_H
#define averagingCurvatureModel_H

#include "divergenceBasedCurvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class averagingCurvatureModel Declaration
\*---------------------------------------------------------------------------*/

class averagingCurvatureModel
:
    public divergenceBasedCurvatureModel
{

private:

    // Private data
    label averagingIterations_;  


    // Private Member Functions
    virtual void computeCurvature(const fvMesh&, const triSurfaceFront&) const;


public:

    TypeName ("averaging");

    // Constructors
    averagingCurvatureModel(const dictionary& configDict);


    //- Destructor
    virtual ~averagingCurvatureModel() = default;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
