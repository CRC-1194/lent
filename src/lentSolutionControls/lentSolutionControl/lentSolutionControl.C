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

\*---------------------------------------------------------------------------*/

#include "lentSolutionControl.H"

#include "surfaceMesh.H"
#include "fvsPatchField.H"

#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(lentSolutionControl, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void lentSolutionControl::checkFluxConvergence()
{
    // No convergence check in the first outer iteration
    if (firstIter())
    {
        volumetricFluxes_.storePrevIter();
        return;
    }

    // Stop flux convergence checks once convergence has been reached
    if (fluxesHaveConverged_)
    {
        return;
    }

    if (corr_ >= maxFluxUpdateIterations_ + 1)
    {
        // Set this flag so convergence of the p-U coupling can be achieved
        fluxesHaveConverged_ = true;
        return;
    }

    // OpenFOAM uses a so-called normalized residual for convergence control
    // which can be comapred in some way with a L1-norm.
    // However, for the fluxes the max-norm is employed
    auto maxRelDelta = max(mag(volumetricFluxes_ - volumetricFluxes_.prevIter())
                        /(max(mag(volumetricFluxes_)).value() + SMALL));
    auto maxAbsDelta = max(mag(volumetricFluxes_ - volumetricFluxes_.prevIter()));

    if (maxRelDelta.value() < relTolVolFlux_ || maxAbsDelta.value() < absTolVolFlux_)
    {
        Info << "Volumetric fluxes have converged:\n\tRelative change: "
             << maxRelDelta.value() << nl
             << "\tAbsolute change: "
             << maxAbsDelta.value()
             << endl;
        fluxesHaveConverged_ = true;
    }
    else
    {
        Info << "\tRelative phi change = " << maxRelDelta.value() << nl
             << "\tAbsolute phi change = " << maxAbsDelta.value()
             << endl;
    }

    volumetricFluxes_.storePrevIter();
}

bool lentSolutionControl::pressureIsConverged() const
{
    return initialPressureResidualBelowThreshold_;
}

void lentSolutionControl::resetControlFlags()
{
    fluxesHaveConverged_ = false;
    updateMomentumFlux_ = true;
    pressureHasConverged_ = false;
    initialPressureResidualBelowThreshold_ = false;
}

void lentSolutionControl::displayConfiguration() const
{
    Info << "\n--------------------------------------------------------" 
         << "\nConfiguration of lent solution procedure parameters"
         << "\n--------------------------------------------------------"
         << endl;

    Info << "\tMaximum number of inner iterations: " << nCorrPISO_
         << "\n\tMaximum number of outer iterations with momentum flux"
         << " update: " << maxFluxUpdateIterations_
         << "\n\tFlux convergence thresholds:"
         << "\n\t\tRelative change: " << relTolVolFlux_
         << "\n\t\tAbsolute change: " << absTolVolFlux_;

    Info << nl << nl;
}

bool lentSolutionControl::evaluatePressureResidual
     (
        const tmp<scalarField> residual,
        const scalarField& rhs
     ) const
{
    scalar minimum = std::numeric_limits<double>::min();
    word resType = dict().subDict("residualControl").subDict("p_rgh").get<word>("resType");
    word norm = dict().subDict("residualControl").subDict("p_rgh").get<word>("norm");
    scalar tolerance = dict().subDict("residualControl").subDict("p_rgh").get<scalar>("tolerance");
    scalarField res{mag(residual.ref())};

    if (resType == "relative")
    {
        res = res/(mag(rhs) + minimum); 
    }

    scalar resNorm = 1.0;

    if (norm == "Linf")
    {
        resNorm = linf(res);
    }
    else if (norm == "L2")
    {
        resNorm = l2(res);
    }
    else if (norm == "L1")
    {
        resNorm = l1(res);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown norm '" << norm << "'.\n"
            << "Use either Linf, L2 or L1."
            << abort(FatalError);
    }

    Info << norm << "(pressure) = " << resNorm << " ("
         << resType << ")" << endl; 

    return resNorm <= tolerance;
}

bool lentSolutionControl::evaluateVelocityResidual
     (
        const tmp<vectorField> residual,
        const vectorField& rhs
     ) const
{
    scalar minimum = std::numeric_limits<double>::min();
    word resType = dict().subDict("residualControl").subDict("U").get<word>("resType");
    word norm = dict().subDict("residualControl").subDict("U").get<word>("norm");
    scalar tolerance = dict().subDict("residualControl").subDict("U").get<scalar>("tolerance");

    vector componentWiseNorm{1,1,1};
    scalar componentMax = 0.0;

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        scalarField resCmpt{mag(residual.ref().component(cmpt))}; 
        scalarField rhsCmpt{rhs.component(cmpt)};

        if (resType == "relative")
        {
            resCmpt = resCmpt/(mag(rhs) + minimum); 
        }

        if (norm == "Linf")
        {
            componentWiseNorm[cmpt] = linf(resCmpt);
        }
        else if (norm == "L2")
        {
            componentWiseNorm[cmpt] = l2(resCmpt);
        }
        else if (norm == "L1")
        {
            componentWiseNorm[cmpt] = l1(resCmpt);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown norm '" << norm << "'.\n"
                << "Use either Linf, L2 or L1."
                << abort(FatalError);
        }

        componentMax = componentWiseNorm[cmpt] > componentMax ?
                         componentWiseNorm[cmpt] : componentMax;
        Info << norm << "(velocity-" << cmpt
             << ") = " << componentWiseNorm[cmpt] << " ("
             << resType << ")" << endl; 
    }

    return componentMax <= tolerance;
}

scalar lentSolutionControl::linf(const scalarField& residual) const
{
    return max(residual);
}

scalar lentSolutionControl::l2(const scalarField& residual) const
{
    return Foam::sqrt(sum(magSqr(residual)/residual.size()));
}

scalar lentSolutionControl::l1(const scalarField& residual) const
{
    return sum(residual)/residual.size();
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
bool lentSolutionControl::read()
{
    pimpleControl::read();

    const dictionary& lentSCDict = dict();

    relTolVolFlux_ = lentSCDict.get<scalar>("phiChangeTolerance");

    // Optional
    if (lentSCDict.found("absPhiChangeTolerance"))
    {
        absTolVolFlux_ = lentSCDict.get<scalar>("absPhiChangeTolerance");
    }

    if (lentSCDict.found("maxFluxUpdateIterations"))
    {
        maxFluxUpdateIterations_ = lentSCDict.get<label>("maxFluxUpdateIterations");
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentSolutionControl::lentSolutionControl(fvMesh& mesh, const surfaceScalarField& phi, const word& dictName)
:
    pimpleControl{mesh, dictName},
    runTime_{mesh.time()},
    volumetricFluxes_{phi},
    relTolVolFlux_{1.0e-8},
    absTolVolFlux_{1.0e-18},
    maxFluxUpdateIterations_{nCorrPIMPLE_},
    fluxesHaveConverged_{false},
    updateMomentumFlux_{true},
    pressureHasConverged_{false},
    initialPressureResidualBelowThreshold_{false}
{
    read();

    displayConfiguration(); 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
label lentSolutionControl::maxFluxUpdateIterations() const
{
    return maxFluxUpdateIterations_;
}

scalar lentSolutionControl::relativeVolumetricFluxTolerance() const
{
    return relTolVolFlux_;
}

scalar lentSolutionControl::absoluteVolumetricFluxTolerance() const
{
    return absTolVolFlux_;
}

bool lentSolutionControl::loop()
{
    read();

    ++corr_;

    label maxOuterIterations = nCorrPIMPLE_;

    // Increase maximum outer iterations in the first time step
    // so the solution algorithm has the chance to correctly initialize the
    // pressure field.
    maxOuterIterations *= runTime_.timeIndex() == 1 ? 10 : 1;

    setFirstIterFlag();

    // Reset control flags in the first iteration of a new time step
    if (firstIter())
    {
        resetControlFlags();
    }

    // Check for convergence or exceeding the maximum of outer iterations
    if (fluxesHaveConverged_ && pressureHasConverged_)
    {
        Info << "LentSC: converged in " << corr_ - 1 << " iterations." << endl;
        corr_ = 0;
        return false;
    }
    else if (corr_ == maxOuterIterations + 1)
    {
        Info << "LentSC: not converged in " << nCorrPIMPLE_ << " iterations" << endl;
        corr_ = 0;
        return false;
    }

    Info << "LentSC: iteration " << corr_ << endl;

    // Checking here for flux convergence ensures a final outer iteration based
    // on the converged fluxes
    checkFluxConvergence();

    return true;
}

bool lentSolutionControl::correctPressure()
{
    bool performInnerIteration = true;

    ++corrPISO_;

    Info << "Pressure correction iteration " << corrPISO_ << endl;

    // Iteration number check
    if (corrPISO_ > nCorrPISO_)
    {
        performInnerIteration = false;
    }
    else if (corrPISO_ > 1)
    {
        // Residual condition
        if (pressureIsConverged())
        {
            performInnerIteration = false;
        }
    }

    if (!performInnerIteration)
    {
        // The condition below is met when the initial residual of the pressure
        // equation is lower than the prescribed threshold from fvSolution
        // in the first corrector iteration
        // ===> pressure has converged
        pressureHasConverged_ = (corrPISO_ == 2);
        corrPISO_ = 0;
    }

    return performInnerIteration;
}

bool lentSolutionControl::explicitVelocityUpdate()
{
    // Explicit update of velocity is only required if the pressure field
    // has changed
    return !pressureIsConverged();
}

bool lentSolutionControl::updateMomentumFlux()
{
    // Final update of momentum fluxes after convergence of the volumetric fluxes
    if (fluxesHaveConverged_)
    {
        updateMomentumFlux_ = false;
        return true;
    }

    return (corr_ <= maxFluxUpdateIterations_ + 1)
                && (firstIter() || updateMomentumFlux_);
}

void lentSolutionControl::solveForVelocity
(
    fvMatrix<vector>& Avelocity
)
{
    bool converged = evaluateVelocityResidual
                     (
                        Avelocity.residual(),
                        Avelocity.source()
                     );

    label solverCalls = 0;
    while (!converged && solverCalls != maxVelocitySolverCalls_)
    {
        Avelocity.solve();

        converged = evaluateVelocityResidual
                    (
                        Avelocity.residual(),
                        Avelocity.source()
                    );
        ++solverCalls;
    }
}

void lentSolutionControl::solveForPressure
(
    fvScalarMatrix& pressureSystem,
    const dictionary& solverDict
)
{
    // Check initial residual
    bool converged = evaluatePressureResidual
                     (
                        pressureSystem.residual(),
                        pressureSystem.source()
                     );

    initialPressureResidualBelowThreshold_ = converged;

    label solverCalls = 0;
    while (!converged && solverCalls != maxPressureSolverCalls_)
    {
        pressureSystem.solve(solverDict);
        
        converged = evaluatePressureResidual
                     (
                        pressureSystem.residual(),
                        pressureSystem.source()
                     );
        ++solverCalls;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
