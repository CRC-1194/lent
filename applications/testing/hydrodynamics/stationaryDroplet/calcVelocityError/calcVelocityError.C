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

Authors
    Tobias Tolle tolle@mathematik.tu-darmstadt.de
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Calculates various norms of the velocity field for each time step.
    Following norms are calculated:
        - one norm
        - two norm
        - maximum norm
    These norms can be calculated for the cell centered velocities as
    well as for the face centered velocites.

    For testcases of the "stationary droplet"-group these norms also 
    represent a measure for the error in velocity.

    See also the paper of Francois et al.:
    (http://dx.doi.org/10.1016/j.jcp.2005.08.004)

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

#include <fstream>
#include <string>

#include "fvCFD.H"
#include "timeSelector.H"

// Functions calculating norms using cell centered velocities
dimensionedScalar calc_one_norm_cc(const volVectorField& U)
{
    dimensionedScalar one_norm = sum(mag(U)) / U.size();

    return one_norm;
}

dimensionedScalar calc_two_norm_cc(const volVectorField& U)
{
    dimensionedScalar two_norm = sqrt(sum(magSqr(U))) / sqrt(dimensionedScalar(U.size()));

    return two_norm;
}

dimensionedScalar calc_maximum_norm_cc(const volVectorField& U)
{
    dimensionedScalar maximum_norm = max(mag(U));

    return maximum_norm;
}

// Function for generating the face centered velocities
void calc_fc_velocities(surfaceScalarField& Uf, surfaceScalarField& phi,
                        const surfaceVectorField& Sf)
{
    Uf = mag(phi) / mag(Sf);
}

// Functions calculating norms using face centered velocities
dimensionedScalar calc_one_norm_fc(const surfaceScalarField& Uf)
{
    dimensionedScalar one_norm = sum(mag(Uf)) / Uf.size();

    return one_norm;
}

dimensionedScalar calc_two_norm_fc(const surfaceScalarField& Uf)
{
    dimensionedScalar two_norm = sqrt(sum(magSqr(Uf))) / sqrt(dimensionedScalar(Uf.size()));

    return two_norm;
}

dimensionedScalar calc_maximum_norm_fc(const surfaceScalarField& Uf)
{
    dimensionedScalar maximum_norm = max(mag(Uf));

    return maximum_norm;
}

// Test if file is empty
bool fileIsEmpty(std::fstream& file)
{
    if (file.tellg() == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Add necessary options and check if they have been set by the user
    argList::addOption
    (
        "errorFile",
        "Path to the error file to which the application appends the results for further analysis."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.optionFound("errorFile"))
    {
        FatalErrorIn("main()")
            << "Please use option '-errorFile' to name the output file."
            << endl << exit(FatalError);
    }

    // Use two separate files for face-centered and cell-centered velocities
    std::string errorFileNameCC = args.optionRead<fileName>("errorFile");
    std::string errorFileNameFC = args.optionRead<fileName>("errorFile");
    errorFileNameCC.append("_cc.dat");
    errorFileNameFC.append("_fc.dat");

    const char* errorFileNameCCPtr = errorFileNameCC.c_str();
    const char* errorFileNameFCPtr = errorFileNameFC.c_str();

    // Open file to write results to
    std::fstream errorFileCC;
    std::fstream errorFileFC;

    std::string header = "# h\t\tdensity ratio\ttime\t\tone-norm\ttwo-norm\tmax-norm";

    errorFileCC.open(errorFileNameCCPtr, std::ios_base::app);
    if (fileIsEmpty(errorFileCC))
    {
        errorFileCC << "# Cell centered velocities\n"
                    << header
                    << std::endl;
    }

    errorFileFC.open(errorFileNameFCPtr, std::ios_base::app);
    if (fileIsEmpty(errorFileFC))
    {
        errorFileFC << "# Face centered velocities\n"
                    << header
                    << std::endl;
    }

    // Set consistent number format
    errorFileCC.precision(4);
    errorFileCC << std::scientific;
    errorFileFC.precision(4);
    errorFileFC << std::scientific;

    // Read varying parameters from dictionaries for unique identification
    // of results
    IOdictionary transportPropertiesDict
    (
        IOobject
        (
            "transportProperties",
            "constant",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    IOdictionary lentSolutionDict
    (
        IOobject
        (
            "lentSolution",
            "system",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary& transportProperties = 
        runTime.lookupObject<dictionary>("transportProperties");
    const dictionary& air = transportProperties.subDict("air");
    const dimensionedScalar rhoAir = air.lookup("rho");
    const dictionary& water = transportProperties.subDict("water");
    const dimensionedScalar rhoWater = water.lookup("rho");

    const scalar densityRatio = rhoWater.value() / rhoAir.value();

    // Get the time directories from the simulation folder using time selector
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    // Initialize the field for face centered velocities
    surfaceScalarField Uf
    (
        IOobject
        (
            "Uf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            dimLength / dimTime,
            0
        )
    );

    // Initialize Sf field
    const surfaceVectorField& Sf = mesh.Sf();

    // Calculate mesh spacing 
    dimensionedScalar h = max((mag(mesh.delta())));

    forAll(timeDirs, timeI)
    {
        // Skip evaluation of initial data
        if (timeI == 0)
        {
            continue;
        }

        // Set the time to the current time directory and update time index
        runTime.setTime(timeDirs[timeI], timeI);

        // Initialize current velocity field from time step directory
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        dimensionedScalar one_norm_cc = calc_one_norm_cc(U);
        dimensionedScalar two_norm_cc = calc_two_norm_cc(U);
        dimensionedScalar maximum_norm_cc = calc_maximum_norm_cc(U);

        // Calculate current face centered velocities und calculate norms
        surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        calc_fc_velocities(Uf, phi, Sf); 

        dimensionedScalar one_norm_fc = calc_one_norm_fc(Uf);
        dimensionedScalar two_norm_fc = calc_two_norm_fc(Uf);
        dimensionedScalar maximum_norm_fc = calc_maximum_norm_fc(Uf);

        if (timeI > 0)
        {
            errorFileCC << h.value() << "\t"
                        << densityRatio << "\t" << runTime.timeName() << "\t\t"
                        << one_norm_cc.value() << "\t"<< two_norm_cc.value()
                        << "\t" << maximum_norm_cc.value() << std::endl;

            errorFileFC << h.value() << "\t"
                        << densityRatio << "\t" << runTime.timeName() << "\t\t"
                        << one_norm_fc.value() << "\t"<< two_norm_fc.value()
                        << "\t" << maximum_norm_fc.value() << std::endl;
        }
    }

    errorFileCC.close();
    errorFileFC.close();
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
