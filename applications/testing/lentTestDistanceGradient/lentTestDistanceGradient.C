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

Application
    testDistanceGradient

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "map.H"
#include "filter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    std::fstream errorFile; 
    errorFile.open("gradientErrors.dat", std::ios_base::app);

    const volScalarField signedDistance(
        IOobject(
            "signedDistance", 
            runTime.timeName(), 
            runTime, 
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const volScalarField alpha1( 
        IOobject(
            "alpha.water", 
            runTime.timeName(), 
            runTime, 
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Cancel out the gradient in the narrow band (interface cells). 
    volScalarField narrowBandFilterField = 
        map(alpha1, 1.0, [](double x) { return  (x > 0) && (x < 1); }); 
    narrowBandFilterField.rename("narrowBandFilterField"); 

    volVectorField distGrad = fvc::grad(signedDistance, "distGrad"); 
    distGrad.writeOpt() = IOobject::AUTO_WRITE; 
    distGrad = distGrad * narrowBandFilterField; 

    volScalarField distGradMag = mag(distGrad); 
    distGradMag.rename("distGradMag"); 
    distGradMag.writeOpt() = IOobject::AUTO_WRITE; 
    distGradMag = distGradMag * narrowBandFilterField;  

    // Filter out the zeros in the bulk. 
    auto distGradMagFilt = filter(distGradMag, [](double x) {return x > 0;}); 

    auto LinfError = mag(1 - distGradMagFilt);  

    dimensionedScalar h = max(mag(mesh.delta())); 
    errorFile << h.value() << " " 
        << max(LinfError)
        << " " << max(distGradMagFilt)
        << " " << min(distGradMagFilt) << "\n";

    runTime.writeNow(); 

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
