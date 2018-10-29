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
bool lentSolutionControl::volFluxesAreConverged()
{
    // Reset the flux convergence flag for a new time step
    if (firstIter())
    {
        fluxesHaveConverged_ = false;
    }
    // To avoid oscillatory behaviour, disable the flux convergence check
    // once the test was successful within a time step
    else if (!fluxesHaveConverged_)
    {
        checkFluxConvergence();
    }

    volumetricFluxes_.storePrevIter();

    return fluxesHaveConverged_;
}

void lentSolutionControl::checkFluxConvergence()
{
    // TODO: check how residuals are computed and mimick this
    // for convergence check of the volumetric fluxes
    // Discuss with Tomislav what would be reasonable. For now, use
    // max norm
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

void lentSolutionControl::displayConfiguration() const
{
    Info << "\n--------------------------------------------------------" 
         << "\nConfiguration of lent solution procedure parameters"
         << "\n--------------------------------------------------------"
         << endl;

    Info << "\tMaximum number of inner iterations: " << nCorrPISO_
         << "\n\tFlux convergence thresholds:"
         << "\n\t\tRelative change: " << relTolVolFlux_
         << "\n\t\tAbsolute change: " << absTolVolFlux_;

    Info << nl << nl;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
bool lentSolutionControl::criteriaSatisfied()
{
    return volFluxesAreConverged() && pimpleControl::criteriaSatisfied();
}

void lentSolutionControl::read()
{
    pimpleControl::read();

    const dictionary& lentSCDict = dict();

    relTolVolFlux_ = lentSCDict.get<scalar>("phiChangeTolerance");
    //absTolVolFlux_ = lentSCDict.get<scalar>("phiChangeTolerance");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentSolutionControl::lentSolutionControl(fvMesh& mesh, const surfaceScalarField& phi, const word& dictName)
:
    pimpleControl{mesh, dictName},
    volumetricFluxes_{phi},
    relTolVolFlux_{1.0e-8},
    absTolVolFlux_{1.0e-18}
{
    read();

    displayConfiguration(); 
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool lentSolutionControl::adaptiveCorrect()
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

bool lentSolutionControl::updateMomentumFlux() const
{
    // Update of momentum fluxes in the first outer iteration
    // has to be ensured
    return firstIter() || !fluxesHaveConverged_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
