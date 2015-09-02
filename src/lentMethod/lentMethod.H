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
    Foam::lentMethod

SourceFiles
    lentMethod.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    A facade class that simplifies the user interface for the LENT hybrid
    front-tracking level-set method by agglomerating a bunch of run-time
    selected small one-algorithm SRP classes.

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


#ifndef lentMethod_H
#define lentMethod_H

#include "regIOobject.H"
#include "typeInfo.H"
#include "fvMesh.H"
#include "triSurfaceFields.H"
#include "lentCommunication.H"
#include "triSurfaceFront.H"
#include "distanceFieldCalculator.H"
#include "frontReconstructionModel.H"
#include "frontReconstructor.H"
#include "triSurfaceFields.H"
#include "frontVelocityCalculator.H"
#include "frontMotionSolver.H"
#include "markerFieldModel.H"
#include "frontSurfaceTensionForceModel.H"

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
    // Method configuration dictionary.
    IOdictionary lentControlDict_;

    // Front <-> volume mesh communication maps. 
    lentCommunication communicationMaps_;   // TODO: Use HashTable for multiphase. TM.

    // Front reconstruction flag.
    mutable bool frontIsReconstructed_;

    // Field calculator
    tmp<distanceFieldCalculator> distanceFieldCalculatorTmp_;

    // Front reconstruction model.
    tmp<frontReconstructionModel> frontReconstructionModelTmp_;

    // Front reconstructor.
    tmp<frontReconstructor> frontReconstructorTmp_;

    // Search algorithm.
    tmp<frontMeshSearch> searchAlgorithmTmp_;

    // Velocity calculator.
    tmp<frontVelocityCalculator> frontVelocityCalculatorTmp_;

    //// Front motion solver.
    tmp<frontMotionSolver> frontMotionSolverTmp_;

    // MarkerField model.
    tmp<markerFieldModel> markerFieldModelTmp_;

    // Surface tension force model.
    tmp<frontSurfaceTensionForceModel> surfaceTensionForceModelTmp_;

    // Member Operators
        void operator=(const lentMethod&); // FIXME: Private in lentCommunication. TM.

public:

    TypeName("lentMethod");

    // Constructors

        lentMethod(
            const triSurfaceMesh& front,
            const fvMesh& mesh,
            word dictName="lentSolution"
        );

        lentMethod(const lentMethod& copy) = default;

    //- Destructor
    virtual ~lentMethod() = default;

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

        void calcMarkerField(volScalarField& markerField) const;

        void reconstructFront(
            triSurfaceFront& front,
            const volScalarField& signedDistance,
            const pointScalarField& pointSignedDistance
        );

        bool isFrontReconstructed() const
        {
            return frontIsReconstructed_;
        }

        void calcFrontVelocity(
            triSurfaceVectorField & frontVelocity,
            const volVectorField& U
        );

        void evolveFront(
            triSurfaceFront& front,
            const triSurfaceVectorField& frontVelocity
        ) const;
        
        tmp<surfaceScalarField> faceSurfaceTensionForce(
            const fvMesh& mesh, 
            const triSurfaceMesh& frontMesh
        ) const
        {
            return surfaceTensionForceModelTmp_->faceSurfaceTensionForce(mesh, frontMesh); 
        }

        tmp<volVectorField> cellSurfaceTensionForce(
            const fvMesh& mesh, 
            const triSurfaceMesh& frontMesh
        ) const
        {
            return surfaceTensionForceModelTmp_->cellSurfaceTensionForce(mesh, frontMesh); 
        }

        const dictionary& dict() const
        {
            return lentControlDict_;
        }

        const frontSurfaceTensionForceModel& surfaceTensionForceModel() const
        {
            return surfaceTensionForceModelTmp_();  
        }

        const frontCurvatureModel& curvatureModel() const
        {
            return surfaceTensionForceModelTmp_->curvatureModel();
        };

        void calcFrontNormals(triSurfaceFront& front) const; 

        // Empty writer: regIOobject is used only for registry queries.
        bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
