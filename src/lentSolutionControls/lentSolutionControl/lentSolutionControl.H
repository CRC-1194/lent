/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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
    Foam::lentSolutionControl

Description
    Class to control the solution procedure of the discretized
    Navier Stokes system for two-phase systems.
    This control class respects and resolves the non-linearity of
    the convective term up to a prescribed tolerance. After reaching this
    tolerance, the momentum fluxes (rhoPhi) are kept constant for the remaining
    outer iterations.
    Update of the momentum fluxes can be controlled by the three parameters
    described below. Both the relative change and the absolute change
    take the maximum norm of the the relative / absolute flux change field.
    The user can also prescribe a maximum number of outer iterations in which
    the momentum fluxes are updated. Beyond this iteration number the
    momentum fluxes are frozen irrespective if the convergence criteria are
    fulfilled.
    Note that the change based criteria are linked with a logical OR. So
    enforcing the fulfillment of both criteria is currently not possible.

    Required parameters:
        - phiChangeTolerance: threshold for realtive flux change convergence

    Optional parameters:
        - absPhiChangeTolerance: threshold for abslute flux change convergence
        - maxFluxUpdateIterations: number of outer iterations after which the
                update of momentum fluxes is stopped.
        
SourceFiles
    lentSolutionControl.C

Authors
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

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

#ifndef lentSolutionControl_H
#define lentSolutionControl_H

#include "fvScalarMatrix.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {


/*---------------------------------------------------------------------------*\
                         Class lentSolutionControl Declaration
\*---------------------------------------------------------------------------*/

class lentSolutionControl
:
    public pimpleControl
{

    // Private data
    
    // Field references
    const Time& runTime_;
    const surfaceScalarField& volumetricFluxes_;

    // User defined flags for solution control
    scalar relTolVolFlux_;
    scalar absTolVolFlux_;
    label maxFluxUpdateIterations_;
    label maxVelocitySolverCalls_ = 10;
    label maxPressureSolverCalls_ = 10;

    // Logical flags which indicate different convergence states
    bool fluxesHaveConverged_;
    bool updateMomentumFlux_;
    bool pressureHasConverged_;
    bool initialPressureResidualBelowThreshold_;

    // Private Member Functions
    void checkFluxConvergence();
    bool pressureIsConverged() const;

    void resetControlFlags();
    void displayConfiguration() const;

    // Residual evaluation for convergence check
    bool evaluatePressureResidual
         (
            const tmp<scalarField> residual,
            const scalarField& rhs
         ) const;
    bool evaluateVelocityResidual
         (
            const tmp<vectorField> residual,
            const vectorField& rhs
         ) const;
    scalar linf(const scalarField& residual) const;
    scalar l2(const scalarField& residual) const;
    scalar l1(const scalarField& residual) const;

    
protected:

    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual bool read() override;


public:

    // Static data members
    
        //- Run-time type information
        TypeName("lentSolutionControl");


    // Constructors
    lentSolutionControl(fvMesh& mesh, const surfaceScalarField& phi, const word& dictName="lentSC");


    //- Destructor
    virtual ~lentSolutionControl() = default;


    // Member Functions

        // Access
        label maxFluxUpdateIterations() const;
        scalar relativeVolumetricFluxTolerance() const;
        scalar absoluteVolumetricFluxTolerance() const;

        // Solution control
        bool loop() override;
        bool correctPressure();
        bool explicitVelocityUpdate();
        bool updateMomentumFlux();

        // Linear solver workaround to allow more control over
        // when the linear systems are considered converged (TT)
        void solveForVelocity
        (
            fvMatrix<vector>& Avelocity
        );
        void solveForPressure
        (
            fvScalarMatrix& pressureSystem,
            const dictionary& solverDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
