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

Author
    Tobias Tolle tolle@csi.tu-darmstadt.de

Description
    Test application for evaluating the quality of the volume fraction
    approximation of a marker field model.

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
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "interfaceProperties.H"
#include "pimpleControl.H"
#include "turbulentTransportModel.H"

#include "analyticalSurface.H"
#include "errorMetrics.H"
#include "lentMarkerFieldTest.H"
#include "lentMethod.H"

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

using namespace FrontTracking;


int main(int argc, char *argv[])
{
    argList::addOption
    (
        "errorFile",
        "Path and name of the file to write the output to"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (!args.optionFound("errorFile"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-errorFile' to set the path and file for "
            << "output." << endl << exit(FatalError);
    }
    
    // Open errorFile in append mode
    const std::string errorFileName = args.optionRead<fileName>("errorFile");
    const char* errorFileNamePtr = errorFileName.c_str();

    std::fstream errorFile;
    errorFile.open(errorFileNamePtr, std::ios_base::app);
    errorFile.precision(4);
    errorFile << std::scientific;

    // Write header
    if (fileIsEmpty(errorFile))
    {
        errorFile << "#grid_spacing\tbounded\tglobal_volume_error\t" 
                  << "interface_volume_error\t"
                  << "local_arithmetic_mean\t"
                  << "local_quadratic_mean\t"
                  << "local_maximum"
                  << std::endl;
    }

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

    tmp<analyticalSurface> analyticalSurfaceTmp( 
        analyticalSurface::New(lent.dict().subDict("analyticalSurface"))
    );

    // Construct front mesh from analytical surface
    frontConstructor frontCon(analyticalSurfaceTmp, mesh);
    front = frontCon.createTriSurface();
    front.write();

    lent.setTriangleCellMapping(frontCon.triangleToCell());

    lent.calcSearchDistances(searchDistanceSqr, pointSearchDistanceSqr);

    // TODO: remove analytical distance calculation once 
    // the distance calculation for boundary mesh points is fixed (TT)
    if (analyticalSurfaceTmp->type() == "Plane")
    {
        frontCon.cellDistance(signedDistance);
        frontCon.pointDistance(pointSignedDistance);
    }
    else
    {
        lent.calcSignedDistances(
            signedDistance,
            pointSignedDistance,
            searchDistanceSqr,
            pointSearchDistanceSqr,
            front
        );
    }

    lent.calcMarkerField(markerField);

    // Write fields for further manual inspection / posprocessing
    markerField.write();
    signedDistance.write();
    pointSignedDistance.write();

    Info << "Start tests...\n" << endl;

    lentMarkerFieldTest test(markerField, front, lent.dict().subDict("markerFieldModel"));

    bool bounded = test.boundedness();
    scalar globalVolumeError = test.globalVolume();
    scalar interfaceVolumeError = test.interfaceVolume();
    List<scalar> localErrors = test.localVolume();

    errorMetrics metrics(localErrors);
    scalar linDev = metrics.arithmeticMeanError();
    scalar quadDev = metrics.quadraticMeanError();
    scalar maxDev = metrics.maximumError();

    Info << "\nTests finished" << endl;

    dimensionedScalar h = max(mag(mesh.delta()));
    errorFile << h.value() << '\t' << bounded << '\t'
              << globalVolumeError << "\t\t"
              << interfaceVolumeError << "\t\t"
              << linDev << "\t\t" << quadDev << "\t\t" << maxDev
              << std::endl;

    errorFile.close();

    Info<< "\nEnd\n" << endl;
    return 0;
};

// ************************************************************************* //
