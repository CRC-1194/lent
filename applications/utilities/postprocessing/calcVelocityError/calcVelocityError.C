/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by original authors
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
    calcVelocityError

Description
    Calculates various norms of the velocity field for each time step.
    Following norms are calculated:
        - one norm
        - two norm
        - maximum norm
    For testcases of the "stationary droplet"-group these norms also 
    represent a measure for the error in velocity.

    TODO:
        * output
        * face centered velocity errors

Authors
    Tobias Tolle tobias.tolle@stud.tu-darmstadt.de

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

// Functions calculating norms using face centered velocities

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Open file to write results to
    std::ofstream output("velocityError.dat", ios_base::out);
    output <<
        "#step\ttime [s]\tone-norm [m/s]\ttwo-norm [m/s]\tmax-norm [m/s]\n";

    // Get the time directories from the simulation folder using time selector
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
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

        dimensionedScalar maximum_norm = calc_maximum_norm_cc(U);
        dimensionedScalar one_norm = calc_one_norm_cc(U);
        dimensionedScalar two_norm = calc_two_norm_cc(U);

//        output << timeI << "\t\t" << runTime.timeName() << "\t\t"
//               << maximum_norm << '\n';
    }

    // Close and write file
    output.close();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
