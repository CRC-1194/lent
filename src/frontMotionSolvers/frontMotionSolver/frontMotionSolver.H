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
    Foam::frontMotionSolver

Description
    Interface for the front motion solution. 

    The motion of the front is defined in a lagrangian way with 

    \partial_t \x_v = \U(\x_v) 

    The solution can be approximated with different (explicity) temporal 
    discretization schemes implemented by classes in the frontMotionSolver
    hierarchy. 

SourceFiles
    frontMotionSolver.C

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef frontMotionSolver_H
#define frontMotionSolver_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"
#include "triSurfaceFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class frontMotionSolver
    :
        public refCount
{

public:

    TypeName ("frontMotionSolver"); 

    declareRunTimeSelectionTable (
        tmp,  
        frontMotionSolver,  
        Dictionary,  
        (
            const dictionary& configDict 
        ),
        (configDict) 
    );

    // Constructors

        frontMotionSolver(const dictionary& configDict);

    // Selectors

        static tmp<frontMotionSolver> New(const dictionary& configDict);

    // Destructor
        virtual ~frontMotionSolver();

    // Member Functions
        
        virtual void evolveFront(
            triSurfaceFront& front,    
            const triSurfaceFrontVectorField& frontVelocity
        ) const = 0;  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
