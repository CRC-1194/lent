/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    setSignedDistnaceFields 

Authors
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

Description
    Pre-processing application that sets two signed distance fields using a 
    surface mesh in the STL format. Used as a pre-processing application for
    the LENT algorithm. Expects to find constant/levelSetFrontDict with 

    narroBandWidth N; 

    defined to define the *square* of the narrow band search radius.

        o cells-to-elements is the signed distance between the surface mesh
          elements and the mesh cell centers

        o points-to-elemetns is the signed distance between the surface mesh
          elements and the mesh points

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "triSurfaceFront.H"
#include "naiveNarrowBandPropagation.H"
#include "TriSurfaceMeshCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace FrontTracking;

// Configure the Method
typedef triSurfaceFront Front;
typedef TriSurfaceMeshCalculator Calculator; 

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Front front(
        IOobject(
            "front.stl",
            "front",
            runTime, 
            IOobject::MUST_READ, 
            IOobject::AUTO_WRITE
        )
    );

    IOdictionary levelSetFrontDict
    (
        IOobject
        (
            "levelSetFrontDict", 
            "constant", 
            runTime, 
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::AUTO_WRITE
        )
    );

    label narrowBandWidth = 
        levelSetFrontDict.lookupOrDefault<label>(
            "narrowBandWidth", 4
        );

    Calculator calc (narrowBandWidth); 

    // Compute the new signed distance field. 
    calc.calcCentresToElementsDistance(
        signedDistance, 
        front,
        naiveNarrowBandPropagation()
    ); 

    signedDistance.write(); 

    calc.calcPointsToElementsDistance(
        pointSignedDistance, 
        front,
        mesh, 
        naiveNarrowBandPropagation()
    ); 

    pointSignedDistance.write(); 

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
