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
    // TODO: check how residuals are computed and mimick this
    // for convergence check of the volumetric fluxes
    // Discuss with Tomislav what would be reasonable. For now, use
    // max norm
    auto maxRelDelta = max(mag(volumetricFluxes_ - volumetricFluxes_.prevIter()))/(max(mag(volumetricFluxes_)).value() + SMALL);

    if (maxRelDelta.value() < relTolVolFlux_)
    {
        Info << "\nVolumetric fluxes have converged. Last relative change: "
             << maxRelDelta.value() << endl;
        return true;
    }
    else
    {
        return false;
    }
}

bool lentSolutionControl::pressureIsConverged()
{
    // Expect no convergence in first inner iteration of first outer iteration
    if (firstIter() && corrPISO_ == 1)
    {
        return false;
    }

    // This is adapted from pimpleControl::criteriaSatisfied()
    const dictionary& solverDict = mesh_.solverPerformanceDict();
    label pressureFieldIndex = -1;
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

void lentSolutionControl::checkParameterCompatibility() const
{
    const word& ut = explicitVelocityUpdateType_;

    if (!(ut == "none" || ut == "adaptive" || ut == "alwaysOn"))
    {
        FatalErrorIn
        (
            "lentSolutionControl::checkParameterCompatibility()"
        )   << "Unknown velocityUpdateType " << ut << "; valid options are"
        << "none, adaptive, alwaysOn" << endl
        << exit(FatalError);
    }

    if (nCorrPISO() > 1 && explicitVelocityUpdateType_ == "none")
    {
        Info << "\nWarning: n = " << nCorrPISO() << " inner iterations set"
             << ", but explicit velocity update is disabled. This makes no sense."
             << endl;
    }

    // Residual controlled velocity solution in combination with explicit velocity
    // update is a bad idea as the explicit update may prevent the fulfillment
    // of the accuracy requirement for the initial residual
    if (!residualControl_.empty())
    {
        for (const auto& residualData : residualControl_)
        {
            if (residualData.name == "U" && explicitVelocityUpdateType_ == "alwaysOn")
            {
                FatalErrorIn
                (
                    "lentSolutionControl::checkParameterCompatibility()"
                )   << "VelocityUpdateType is " << ut << "; this will not work with"
                << " a residual controlled solution of the velocity." << endl
                << exit(FatalError);
            }
        }
    }
}

void lentSolutionControl::displayConfiguration() const
{
    Info << "\n--------------------------------------------------" 
         << "\nConfiguration of solution procedure"
         << "\n--------------------------------------------------"
         << endl;

    Info << "PIMPLE configuration:\n";

    if (residualControl_.empty())
    {
        Info<< ": no residual control data found. "
            << "Calculations will employ " << nCorrPIMPLE_
            << " outer corrector loops" << nl;
    }
    else
    {
        Info<< ": max outer iterations = " << nCorrPIMPLE_ << nl;

        for (const fieldData& ctrl : residualControl_)
        {
            Info<< "    field " << ctrl.name << token::TAB
                << ": relTol " << ctrl.relTol
                << ", tolerance " << ctrl.absTol
                << nl;
        }
    }

    Info << "LentSolutionControl specific parameters:\n"
         << "\tResolve convective non-linearity: " << resolveConvectiveNonlinearity_
         << "\n\tUse adaptive inner iteration strategy: " << adaptiveInnerIteration_
         << "\n\tExplicit velocity update: " << explicitVelocityUpdateType_;

    if (explicitVelocityUpdateType_ == "adaptive")
    {
        Info << "\n\tRelative tolerance for volumetric flux change: "
             << relTolVolFlux_
             << "\nExplicit velocity update will be disabled when this"
             << " tolerance is reached.\n";
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void lentSolutionControl::read()
{
    pimpleControl::read();

    const dictionary& lentSCDict = dict();

    resolveConvectiveNonlinearity_ = lentSCDict.lookupType<bool>("updateMomentumFlux");
    adaptiveInnerIteration_ = lentSCDict.lookupType<bool>("adaptiveInnerIteration");
    explicitVelocityUpdateType_ = lentSCDict.lookupType<word>("velocityUpdateType");

    if (explicitVelocityUpdateType_ == "adaptive")
    {
        relTolVolFlux_ = readScalar(lentSCDict.lookup("phiChangeTolerance"));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentSolutionControl::lentSolutionControl(fvMesh& mesh, const surfaceScalarField& phi, const word& dictName)
:
    pimpleControl{mesh, dictName},
    volumetricFluxes_{phi},
    resolveConvectiveNonlinearity_{true},
    adaptiveInnerIteration_{true},
    explicitVelocityUpdateType_{"adaptive"},
    relTolVolFlux_{1.0e-8},
    explicitVelocityUpdate_{true}
{
    read();

    // TODO: add some checks for contradicting parameter sets.
    // E.g. it makes no sense to have more than one inner iteration without
    // the explicit velocity update.
    checkParameterCompatibility();

    // Catch those combinations here.
    // TODO: output information of solution control configuration
    // analogue to pimpleControl
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

    // Check adaptive conditions
    if (adaptiveInnerIteration_)
    {
        // Residual condition
        if (pressureIsConverged())
        {
            performInnerIteration = false;
        }
        // Without explicit velocity update, further inner iterations have no
        // effect
        else if (!explicitVelocityUpdate_ && corrPISO_ > 1)
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
    if (explicitVelocityUpdateType_ == "none")
    {
        explicitVelocityUpdate_ = false;
    }
    else if (explicitVelocityUpdateType_ == "alwaysOn")
    {
        explicitVelocityUpdate_ = true;
    }
    else if (explicitVelocityUpdateType_ == "adaptive")
    {
        // Explicit update must be re-enabled for each new time step
        if (firstIter())
        {
            explicitVelocityUpdate_ = true;
        }
        else
        {
            if (explicitVelocityUpdate_)
            {
                explicitVelocityUpdate_ = !volFluxesAreConverged();
            }
        }

        volumetricFluxes_.storePrevIter();
    }

    return explicitVelocityUpdate_;
}

bool lentSolutionControl::updateMomentumFlux() const
{
    // Update of momentum fluxes in the first outer iteration
    // has to be ensured
    return firstIter() || resolveConvectiveNonlinearity_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
