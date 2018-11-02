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

    // OpenFOAM uses a so-called normalized residual for convergence control
    // which can be comapred in some way with a L1-norm.
    // However, for the fluxes the max-norm is employed
    auto maxRelDelta = max(mag(volumetricFluxes_ - volumetricFluxes_.prevIter())
                        /(max(mag(volumetricFluxes_)).value() + SMALL));
    auto maxAbsDelta = max(mag(volumetricFluxes_ - volumetricFluxes_.prevIter()));

    if (maxRelDelta.value() < relTolVolFlux_ || maxAbsDelta.value() < absTolVolFlux_)
    {
        Info << "\nVolumetric fluxes have converged:\n\tLast relative change: "
             << maxRelDelta.value() << nl
             << "\tLast absolute change: "
             << maxAbsDelta.value() << nl
             << endl;
        fluxesHaveConverged_ = true;
    }
    else
    {
        Info << "\nPhi change = " << maxRelDelta.value() << nl
             << "max(mag(phi)) = " << max(mag(volumetricFluxes_)).value() << nl
             << "Abs phi change = " << maxAbsDelta.value() << nl
             << endl;
    }

    volumetricFluxes_.storePrevIter();
}

bool lentSolutionControl::pressureIsConverged() const
{
    // Expect no convergence in first inner iteration of first outer iteration
    if (firstIter() && corrPISO_ == 1)
    {
        return false;
    }

    // This is adapted from pimpleControl::criteriaSatisfied()
    const dictionary& solverDict = mesh_.solverPerformanceDict();
    label pressureFieldIndex = -2;
    Pair<scalar> residuals{GREAT, GREAT};

    forAllConstIters(solverDict, iter)
    {
        const entry& solverPerfDictEntry = *iter;

        const word& fieldName = solverPerfDictEntry.keyword();

        // TODO: hard code name of pressure field? (TT)
        if (fieldName == "p_rgh")
        {
            pressureFieldIndex = applyToField(fieldName);
            residuals = maxResidual(solverPerfDictEntry);
            break;
        }
    }

    if (residuals.last() < residualControl_[pressureFieldIndex].absTol)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void lentSolutionControl::resetControlFlags()
{
    fluxesHaveConverged_ = false;
    updateMomentumFlux_ = true;
    pressureHasConverged_ = false;
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void lentSolutionControl::read()
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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentSolutionControl::lentSolutionControl(fvMesh& mesh, const surfaceScalarField& phi, const word& dictName)
:
    pimpleControl{mesh, dictName},
    volumetricFluxes_{phi},
    relTolVolFlux_{1.0e-8},
    absTolVolFlux_{1.0e-18},
    maxFluxUpdateIterations_{nCorrPIMPLE_},
    fluxesHaveConverged_{false},
    updateMomentumFlux_{true},
    pressureHasConverged_{false}
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
    else if (corr_ == nCorrPIMPLE_ + 1)
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
