/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::frontReconstructionModel

Description
    Interface for the reconstruction test: LENT reconstruction is to be performed
    only when topological changes are happening and / or the quality of the front
    has decreased. How and when exactly is what this class hierarchy decides on.

SourceFiles
    frontReconstructionModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef frontReconstructionModel_H
#define frontReconstructionModel_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontReconstructionModel Declaration
\*---------------------------------------------------------------------------*/

class frontReconstructionModel
    :
        public refCount
{

public:

    TypeName ("frontReconstructionModel"); 

    declareRunTimeSelectionTable (
        tmp,  
        frontReconstructionModel,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors

        explicit frontReconstructionModel(const dictionary& configDict);

    // Selectors

        static tmp<frontReconstructionModel> New(const dictionary& configDict);

    // Destructor
        virtual ~frontReconstructionModel();

    // Member Functions
        
        virtual bool reconstructionRequired(
            const triSurfaceFront& front, 
            const volScalarField& signedDistance 
        ) const = 0;  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
