/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    testDistanceGradient

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    std::fstream errorFile; 
    errorFile.open("gradientErrors.dat", std::ios_base::app);

    volScalarField signedDistance(
        IOobject(
            "signedDistance", 
            runTime.timeName(), 
            runTime, 
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField alpha1( 
        IOobject(
            "alpha.water", 
            runTime.timeName(), 
            runTime, 
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField distGrad = fvc::grad(signedDistance, "distanceGradient"); 

    volScalarField distGradMag = mag(distGrad); 

    distGradMag *= map(alpha1, 1.0, [](double x) { return  (x > 0) && (x < 1); }); 

    Info << max(distGradMag).value() << " " << min(distGradMag).value() << endl;

    errorFile  << max(distGradMag).value() << " " << min(distGradMag).value() << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
