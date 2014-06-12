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
    Foam::eulerFrontMotionSolver

Description
    Evolve the front using the first order accurate Euler temporal discretization
    scheme. 

SourceFiles
    eulerFrontMotionSolver.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#ifndef eulerFrontMotionSolver_H
#define eulerFrontMotionSolver_H

#include "frontMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class eulerFrontMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class eulerFrontMotionSolver
    :
        public frontMotionSolver 
{

public:

    TypeName ("Euler"); 

    // Constructors

        eulerFrontMotionSolver(const dictionary& configDict);

    // Destructor
        virtual ~eulerFrontMotionSolver();

    // Member Functions
        
        void evolveFront(
            triSurfaceFront& front, 
            const triSurfaceFrontVectorField& frontVelocity 
        ) const;  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
