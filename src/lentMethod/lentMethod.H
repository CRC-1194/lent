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
    Foam::lentMethod

SourceFiles
    lentMethod.C

Authors
    Tomislav Maric (maric@csi.tu-darmstadt.de)
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    A facade class that simplifies the user interface for the LENT hybrid
    front-tracking level-set method by agglomerating a bunch of run-time
    selected small one-algorithm SRP classes.

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


#ifndef lentMethod_H
#define lentMethod_H

#include "regIOobject.H"
#include "typeInfo.H"
#include "fvMesh.H"
#include "lentCommunication.H"
#include "triSurfaceFront.H"
#include "distanceFieldCalculator.H"
#include "frontReconstructionModel.H"
#include "frontReconstructor.H"
#include "triSurfaceFields.H"
#include "frontMotionSolver.H"
#include "markerFieldModel.H"
#include "frontSurfaceTensionForceModel.H"
#include "lentInterpolation.H"
#include "triSurfaceFrontFields.H"
// TODO: move this to a place more suitable, e.g. frontReconstructionModel
#include "reconstructionHistory.H"

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
    lentCommunication communicationMaps_;   

    // Front reconstruction flag.
    mutable bool frontIsReconstructed_;

    // Field calculator
    tmp<distanceFieldCalculator> distanceFieldCalculatorTmp_;

    // Front reconstruction model.
    tmp<frontReconstructionModel> frontReconstructionModelTmp_;

    // Front reconstructor.
    tmp<frontReconstructor> frontReconstructorTmp_;

    // Front motion solver.
    tmp<frontMotionSolver> frontMotionSolverTmp_;

    // MarkerField model.
    tmp<markerFieldModel> markerFieldModelTmp_;

    // Surface tension force model.
    tmp<frontSurfaceTensionForceModel> surfaceTensionForceModelTmp_;
    
    // Algorithm collection for front smoothing
    frontSmoother frontSmoother_;

    // Recorder for front reconstructions and smoothing
    reconstructionHistory reconstructionHistory_;

    // Hack to enable smoothing without reconstruction
    Switch useReconstruction_;

    // Member Operators
        void operator=(const lentMethod&); // FIXME: Private in lentCommunication. TM.

public:

    TypeName("lentMethod");

    // Constructors

        lentMethod(
            const triSurfaceFront& front,
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
            triSurfaceFrontPointVectorField & frontVelocity,
            const volVectorField& U
        );

        void evolveFront(
            triSurfaceFront& front,
            const volVectorField& cellVelocity
        );
        
        tmp<surfaceScalarField> faceSurfaceTensionForce(
            const fvMesh& mesh, 
            const triSurfaceFront& frontMesh
        ) const
        {
            return surfaceTensionForceModelTmp_->faceSurfaceTensionForce(mesh, frontMesh); 
        }

        tmp<volVectorField> cellSurfaceTensionForce(
            const fvMesh& mesh, 
            const triSurfaceFront& frontMesh
        ) const
        {
            return surfaceTensionForceModelTmp_->cellSurfaceTensionForce(mesh, frontMesh); 
        }

        tmp<fvMatrix<vector>> surfaceTensionImplicitPart(
            const volVectorField& velocity,
            const volScalarField& markerField,
            const triSurfaceFront& front
        ) const
        {
            return surfaceTensionForceModelTmp_->surfaceTensionImplicitPart(velocity, markerField, front);
        }

        const dictionary& dict() const
        {
            return lentControlDict_;
        }

        const frontSurfaceTensionForceModel& surfaceTensionForceModel() const
        {
            return surfaceTensionForceModelTmp_();  
        }

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
