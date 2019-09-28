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

    adjustMaxCoTimeStep 

Description
    
    Adjusts the initial time step for a simulation case based on computed
    maximal Courant number.

Authors
    Tomislav Maric tomislav@sourceflux.de 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// Advection test field models.
#include "divFreeFieldModel.H"

scalar maxCoNum(const auto& phi)
{
    const auto& mesh = phi.mesh(); 
    const auto& runTime = mesh.time(); 

    #include "lentCourantNo.H"

    Info << "dictionary maxCo = " << maxCo << endl;

    return CoNum;  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:


int main(int argc, char *argv[])
{
    #include "setRootCase.H" 
    #include "createTime.H" 
    #include "createMesh.H" 

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    volVectorField U 
    (
        IOobject
        (
            "U", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
        mesh, 
        dimensionedVector("U", dimLength / dimTime, vector(0,0,0))
    ); 
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
        mesh, 
        dimensionedScalar("phi", dimVolume / dimTime, 0)
    ); 

    IOdictionary controlDict 
    (
        IOobject
        (
            "controlDict", 
            runTime.system(), 
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    const scalar CoDict = readScalar(controlDict.lookup("maxCo")); 
    const scalar deltaTdict = readScalar(controlDict.lookup("deltaT")); 

    // Run function objects in case they set the volumetric flux field.
    // Select the divergence free velocity/flux model. 
    autoPtr<divFreeFieldModel> divFree = divFreeFieldModel::New(runTime);
    divFree->execute();  

    auto CoNum = maxCoNum(phi);  

    Info << "CoNum = " << CoNum << endl;
    if (CoNum > SMALL)
    {
        // Scale the initial time step
        double deltaTdictNew = deltaTdict * (CoDict / CoNum);

        Info << "New time step = " << deltaTdictNew << endl;

        unsigned int nTimeSteps = floor(runTime.endTime().value() / deltaTdictNew);

        if (nTimeSteps % 2 != 0)
            nTimeSteps += 1; 

        Info << "nTimesteps = " << nTimeSteps << endl;

        deltaTdictNew = runTime.endTime().value() / nTimeSteps; 

        Info << "Rounded time step = " << deltaTdictNew << endl; 

        runTime.setDeltaT(deltaTdictNew); 

        // Check the new Co and deltaT. 
        maxCoNum(phi); 

        controlDict.add("deltaT", deltaTdictNew, true); 
        const regIOobject& reg = controlDict;
        reg.write(); 

        Info<< "\nEnd\n" << endl;
    }

    return 0;
}


// ************************************************************************* //
