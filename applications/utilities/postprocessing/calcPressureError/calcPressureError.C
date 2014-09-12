/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x
    \\  /    A nd           | Copyright held by original author
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

Application
    calcPressureError

Description
    Calculate three pressure errors as defined in the paper of Francois
    (see http://dx.doi.org/10.1016/j.jcp.2005.08.004) for every time step

    TODO:
        * Way to distinguish between circle and sphere without user input
          --> based on presence of empty patches?
          For now an user argument is used

Author
    Tobias Tolle tobias.tolle@stud.tu-darmstadt.de

\*---------------------------------------------------------------------------*/
#include <fstream>
#include <string>

#include "fvCFD.H"
#include "timeSelector.H"

// Calculate difference of mean pressure inside and outside of the drop
scalar calc_deltaP_total(volScalarField& P, const volVectorField& cellCenters,
                         const scalar& radius, const vector& center)
{
    scalar deltaP_total = 0;

    scalar P_in = 0;
    scalar P_out = 0;
    scalar n_drop_cells = 0;
    vector distance (0, 0, 0);

    // Sum up pressure for averaging, distinguish between drop and
    // background fluid
    forAll(P, I)
    {
        distance = center - cellCenters[I];

        if (mag(distance) > radius)
        {
            P_out += P[I];
        }
        else
        {
            P_in += P[I];
            n_drop_cells++;
        }
    }

    forAll(P.boundaryField(), I)
    {
        scalarField& PBoundaryField = P.boundaryField()[I];

        forAll(PBoundaryField, J)
        {
            P_out += PBoundaryField[J];
        }
    }

    P_in /= n_drop_cells;
    P_out = P_out / (P.size() - n_drop_cells);

    deltaP_total = P_in - P_out;

    return deltaP_total;
}

// Calculate difference of mean pressure inside and outside of the drop
// without the transition zone around the interface
scalar calc_deltaP_partial(volScalarField& P, const volVectorField& cellCenters,
                           const scalar& radius, const vector& center)
{
    scalar deltaP_partial = 0;

    scalar P_in = 0;
    scalar P_out = 0;
    scalar n_drop_cells = 0;
    scalar n_trans_cells = 0;
    vector distance (0, 0, 0);

    // Sum up pressure for averaging, distinguish between drop and background
    // fluid. Cells in the transition region are not considered in this case
    forAll(P, I)
    {
        distance = center - cellCenters[I];

        if (mag(distance) <= 0.5*radius)
        {
            P_in += P[I];
            n_drop_cells++;
        }
        else if (mag(distance) >= 1.5*radius)
        {
            P_out += P[I];
        }
        else
        {
            n_trans_cells++;
        }
    }

    forAll(P.boundaryField(), I)
    {
        scalarField& PBoundaryField = P.boundaryField()[I];

        forAll(PBoundaryField, J)
        {
            P_out += PBoundaryField[J];
        }
    }

    P_in /= n_drop_cells;
    P_out = P_out / (P.size() - n_drop_cells - n_trans_cells);

    deltaP_partial = P_in - P_out;

    return deltaP_partial;
}

// Calculate maximum pressure difference
scalar calc_deltaP_max(volScalarField& P)
{
    dimensionedScalar deltaP_max = max(P) - min(P);

    return deltaP_max.value();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Add necessary options and check if they have been set by the user
    argList::addOption
    (
        "errorFile",
        "Path and name of the file to write the output to"
    );

    argList::addOption
    (
        "radius",
        "Radius of the circle/sphere used in your testcase"
    );

    argList::addOption
    (
        "center",
        "Vector to the center of the circle/sphere"
    );

    argList::addOption
    (
        "shape",
        "Specify if circle or sphere"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.optionFound("errorFile"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-errorFile' to set the path and file for "
            << "output." << endl << exit(FatalError);
    }

    if (!args.optionFound("radius"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-radius' to specify the radius of your "
            << "testcase." << endl << exit(FatalError);
    }

    if (!args.optionFound("center"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-center' to specify the circle/sphere "
            << "center." << endl << exit(FatalError);
    }

    if (!args.optionFound("shape"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Pleasue use option '-shape' to specify if calculation is done"
            << " for a circle or a sphere" << endl << exit(FatalError);
    }

    // Open errorFile in append mode
    const std::string errorFileName = args.optionRead<fileName>("errorFile");
    const char* errorFileNamePtr = errorFileName.c_str();

    std::fstream errorFile;
    errorFile.open(errorFileNamePtr, std::ios_base::app);
    errorFile.precision(4);

    // Write header
    errorFile << "# step\ttime[s]\terror p_total\terror p_partial\t"
              << "error p_max\n"; 

    const scalar radius = args.optionRead<scalar>("radius");
    const vector center = args.optionRead<vector>("center");
    const word shape = args.optionRead<word>("shape");
    const volVectorField& C = mesh.C();
    
    // Read surface tension coefficient to calculate exact pressure jump
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
    const dimensionedScalar sigma(transportPropertiesDict.lookup("sigma"));
    scalar deltaP_exact = 0;

    // Set pressure jump according to selected shape
    if (shape == "circle")
    {
        deltaP_exact = sigma.value() / radius;
    }
    else if (shape == "sphere")
    {
        deltaP_exact = sigma.value() * 2 / radius;
    }
    else
    {
        Info << "Shape " << shape << "unknown, assuming circle\n" << endl;
        deltaP_exact = sigma.value() / radius;
    }

    // Get the time directories from the simulation folder using time selector
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        // Read pressure field of current time step
        volScalarField P
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        // Calculate numeric pressure differences and errors
        scalar deltaP_total = calc_deltaP_total(P, C, radius, center);
        scalar deltaP_partial = calc_deltaP_partial(P, C, radius, center);
        scalar deltaP_max = calc_deltaP_max(P);

        scalar error_total = mag(mag(deltaP_total) - mag(deltaP_exact)) / mag(deltaP_exact);
        scalar error_partial = mag(mag(deltaP_partial) - mag(deltaP_exact)) / mag(deltaP_exact);
        scalar error_max = mag(mag(deltaP_max) - mag(deltaP_exact)) / mag(deltaP_exact);

        // Write errors to file, ignore initial condition
        if (timeI > 0)
        {
            errorFile << timeI << "\t\t" << runTime.timeName() << "\t"
                      << error_total << "\t\t\t" << error_partial << "\t\t\t"
                      << error_max << "\n";
        }

    }

    // Close errorFile
    errorFile.close();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
