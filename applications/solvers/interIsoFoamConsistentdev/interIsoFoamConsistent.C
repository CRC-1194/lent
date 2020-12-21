/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016 DHI
    Copyright (C) 2017 OpenCFD Ltd.
    Copyright (C) 2018 Johan Roenby
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

Application
    interIsoFoam

Group
    grpMultiphaseSolvers

Description
    Solver derived from interFoam for two incompressible, isothermal immiscible
    fluids using the isoAdvector phase-fraction based interface capturing
    approach, with optional mesh motion and mesh topology changes including
    adaptive re-meshing.

    Reference:
    \verbatim
        Roenby, J., Bredmose, H. and Jasak, H. (2016).
        A computational method for sharp interface advection
        Royal Society Open Science, 3
        doi 10.1098/rsos.160405
    \endverbatim

    isoAdvector code supplied by Johan Roenby, STROMNING (2018)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "isoCutFace.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using isoAdvector phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing.\n"
        "The solver is derived from interFoam"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();

    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors) // in fvSolution dictory keyword moveMeshOuterCorrectors, used to ??
            {
                mesh.update();

                if (mesh.changing())
                {

                    gh = (g & mesh.C()) - ghRef; // mesh.C-->mesh center position + vector
                    ghf = (g & mesh.Cf()) - ghRef; // mesh.Cf -->mesh face center position + vector

                    MRF.update();

                    if (correctPhi)  // when the mesh change, correct flux
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        mixture.correct();
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // read alpha1 name and number of alpha subcycles
            #include "alphaControls.H"

            // it is beneficial to solve the phase fraction equation in several subcycles within a single time step.
            //DOI:https://doi.org/10.1103/PhysRevE.79.036306
            #include "alphaEqnSubCycle.H" 

            // mixture--> an object of immiscibleIncompressibleTwoPhaseMixture class, which inherits from 
            // 1.incompressibleTwoPhaseMixture-->correct laminar viscosity to keep boundness of alpha
            // 2.interfaceProperties--> correct the alpha boundary condition for dynamic contact angle.  
            mixture.correct();
//------------------------------added contents------------------------------------//             
            #include "rhofs.H"
            scalarField pfvalue = volPointInterpolation::New(mesh).interpolate(alpha1);
            isoCutFace isoface(mesh, pfvalue);
            dimensionedScalar areaDim("areaDim",dimensionSet(0,2,0,0,0,0,0),1.0);
            surfaceScalarField alphaFace = mag(isoface.subFaceArea())*areaDim/mesh.magSf();
            surfaceScalarField rhof = alphaFace*rho1f + (1.0 - alphaFace)*rho2f;
            surfaceScalarField muf = alphaFace*rho1f*nu1 + (1.0 - alphaFace)*rho2f*nu2;
            surfaceScalarField rhoPhi = rhof * phi;                        

            fvScalarMatrix rhoEqn
            (
                 fvm::ddt(rho) + fvc::div(rhoPhi)
             );
            rhoEqn.solve();
//----------------------------------End-------------------------------------------//
            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
