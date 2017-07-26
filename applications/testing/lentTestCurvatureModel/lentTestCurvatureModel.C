/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by orignal authors.
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
    lentTestCurvatureModels 

Description
    Test LENT curvature models against the exact curvature. 

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "lentMethod.H"
#include "volMesh.H"
#include "map.H"
#include <fstream>

using namespace FrontTracking;

// Test if file is empty: for an empty file the current position in a stream in
// append-mode is 0
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
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    std::fstream errorFile; 

    errorFile.open("curvatureErrors.dat", std::ios_base::app);

    triSurfaceFront front(
        IOobject(
            "front",
            "front",
            runTime,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    lentMethod lent(front, mesh);

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    runTime++;

    lent.reconstructFront(front, signedDistance, pointSignedDistance);

    lent.calcSignedDistances(
        signedDistance,
        pointSignedDistance,
        searchDistanceSqr,
        pointSearchDistanceSqr,
        front
    );

    const dictionary& lentDict = lent.dict(); 

    tmp<frontCurvatureModel> exactCurvatureModelTmp = frontCurvatureModel::New(lentDict.subDict("exactCurvatureModel")); 
    const frontCurvatureModel& exactCurvatureModel = exactCurvatureModelTmp.ref(); 
    tmp<volScalarField> cellCurvatureExactTmp = exactCurvatureModel.cellCurvature(mesh,front);  
    volScalarField& exactCurvature = cellCurvatureExactTmp.ref();  

    const auto& surfaceTensionDict = lentDict.subDict("surfaceTensionForceModel");
    tmp<frontCurvatureModel> numericalCurvatureModelTmp = frontCurvatureModel::New(surfaceTensionDict.subDict("curvatureModel"));  
    const frontCurvatureModel& numericalCurvatureModel = numericalCurvatureModelTmp.ref(); 
    tmp<volScalarField> numericalCurvatureTmp = numericalCurvatureModel.cellCurvature(mesh,front);  
    volScalarField& numericalCurvature = numericalCurvatureTmp.ref();  
    numericalCurvature.rename("numericalCurvature"); 
    numericalCurvature.writeOpt() = IOobject::AUTO_WRITE; 

    // Use the gradient of the markerField as filter since this gives the cells
    // where the curvature is used when the CSF model is used for surface tension
    dimensionedScalar dSMALL("SMALL", pow(dimLength,-1), SMALL);
    volScalarField onesFilter = pos(mag(fvc::grad(markerField)) - 1e4*dSMALL);
    onesFilter.rename("ones"); 
    onesFilter.writeOpt() = IOobject::AUTO_WRITE; 
    numericalCurvature *= onesFilter; 
    exactCurvature *= onesFilter; 

    Info << "\nMaximum numerical curvature: " << max(numericalCurvature).value()
         << "\nMinimum numerical curvature: " << min(numericalCurvature).value() << endl; 

    volScalarField LinfField ("LinfCurvatureErr", mag(exactCurvature - numericalCurvature)/(dSMALL + mag(exactCurvature))); 
    LinfField.writeOpt() = IOobject::AUTO_WRITE; 
    dimensionedScalar Linf = max(LinfField);

    // Write results
    if (fileIsEmpty(errorFile))
    {
        errorFile << "# mesh spacing | Linf relative | exact model | numerical model\n";
    }

    dimensionedScalar h = max(mag(mesh.delta())); 
    errorFile << h.value() << " " << Linf.value() << " " << exactCurvatureModel.type()
              << " " << numericalCurvatureModel.type() << std::endl;

    front.write();
    runTime.writeNow();

    Info <<"Maximal relative curvature error = " << Linf.value() << '\n' << endl;

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
