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
    Foam::lentMethod 

Description
    A facade class that simplifies the user interface for the LENT hybrid 
    front-tracking level-set method by agglomerating a bunch of run-time 
    selected small one-algorithm SRP classes. 

SourceFiles
    lentMethod.C

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#ifndef lentMethod_H
#define lentMethod_H

#include "regIOobject.H"
#include "typeInfo.H"
#include "triSurfaceFront.H"
#include "fvMesh.H"

#include "lentDistanceFieldCalculator.H"
#include "frontReconstructionModel.H"
#include "frontReconstructor.H"
#include "triSurfaceFrontFields.H"
#include "frontVelocityCalculator.H"
//#include "frontMotionSolver.H"
#include "heavisideModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

/*---------------------------------------------------------------------------*\
                         Class lentMethod Declaration
\*---------------------------------------------------------------------------*/


class lentMethod
:
    public regIOobject
{
    // Triangles-cells map.
    DynamicList<label> meshCells_; 

    // Lent method configuration dictionary.
    IOdictionary lentControlDict_;  

    // Field calculator
    tmp<lentDistanceFieldCalculator> lentDistanceFieldCalculatorTmp_; 

    // Front reconstruction model. 
    tmp<frontReconstructionModel> frontReconstructionModelTmp_;  

    // Front reconstructor.
    tmp<frontReconstructor> frontReconstructorTmp_; 

    // Velocity calculator. 
    tmp<frontVelocityCalculator> frontVelocityCalculatorTmp_; 

    //// Front motion solver. 
    //tmp<frontMotionSolver> frontMotionSolver_; 

    // Heaviside model. 
    tmp<heavisideModel> heavisideModelTmp_; 
    
public:

    TypeName("lentMethod"); 

    // Constructors

        lentMethod(
            const triSurfaceFront& front, 
            const fvMesh& mesh, 
            word dictName="lentSolution"
        );

        lentMethod(const lentMethod& copy);

    //- Destructor
    virtual ~lentMethod();

    // Member Functions
    
        // Calculations
        virtual void calcSearchDistances(
            volScalarField& searchDistanceSqr, 
            pointScalarField& pointSearchDistanceSqr
        ); 

        void calcSignedDistances(
            volScalarField& signedDistance, 
            pointScalarField& pointSignedDistance,
            const volScalarField& searchDistanceSqr, 
            const pointScalarField& pointSearchDistanceSqr,
            const triSurfaceFront& front
        );

        void calcHeaviside(
            volScalarField& heaviside,
            const volScalarField& signedDistance,
            const volScalarField& searchDistanceSqr
        ) const; 

        void reconstructFront(
            triSurfaceFront& front,
            const volScalarField& signedDistance,
            const pointScalarField& pointSignedDistance
        ); 

        void calcFrontVelocity(
            triSurfaceFrontVectorField& frontVelocity, 
            const volVectorField& U 
        ); 

        // Empty writer: regIOobject is used only for registry queries. 
        bool writeData(Ostream& os) const; 

    // Member Operators
    //
        void operator=(const lentMethod&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
