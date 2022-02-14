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
    Foam::frontMotionSolver

SourceFiles
    frontMotionSolver.C

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Interface for the Front Tracking advection. 

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


#ifndef frontMotionSolver_H
#define frontMotionSolver_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"
#include "triSurfaceFrontFields.H"
#include "lentInterpolation.H"
#include "frontSmoother.H"

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
    tmp<volVectorField> cellDisplacementTmp_; 

    tmp<triSurfaceFrontPointVectorField> frontDisplacementTmp_; 

    lentInterpolation interpolation_;  


protected:

    triSurfaceFrontPointVectorField& frontDisplacements() 
    {
        return frontDisplacementTmp_.ref();
    }

    const triSurfaceFrontPointVectorField& frontDisplacements() const 
    {
        return frontDisplacementTmp_();
    }

    volVectorField& cellDisplacements() 
    {
        return cellDisplacementTmp_.ref();
    }

    const volVectorField& cellDisplacements() const
    {
        return cellDisplacementTmp_();
    }

    void initDisplacements(
        const triSurfaceFront& front, 
        const volVectorField& cellVelocity
    ); 

    virtual void calcCellDisplacement(const volVectorField& cellVelocity) = 0;  

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
    )

    // Constructors
        explicit frontMotionSolver(const dictionary& configDict);

    // Selectors
        static tmp<frontMotionSolver> New(const dictionary& configDict);

    // Destructor
        virtual ~frontMotionSolver() = default;

    // Member Functions
        virtual void evolveFront(
            triSurfaceFront& front,
            const volVectorField& cellVelocity 
        ) final; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
