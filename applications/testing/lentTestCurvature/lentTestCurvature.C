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
    calcExactCurvature 

Description
    Calculates exact curvature for a circle. 

    TODO:  
        * sphere, ellipsoid
        * run-time selection of shapes

Authors
    Tobias Tolle tobias.tolle@stud.tu-darmstadt.de  
    Tomislav Maric maric@csi.tu-darmstadt.de

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

scalar circle_curvature(const point& p, const point& circleCenter)
{
    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    surfaceScalarField faceCurvature  
    (
        IOobject
        (
            "faceCurvature",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh, 
        dimensionedScalar(
            "zero", 
            pow(dimLength, -1), 
            0
        ) 
    );

    volScalarField cellCurvature  
    (
        IOobject
        (
            "cellCurvature",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh, 
        dimensionedScalar(
            "zero", 
            pow(dimLength, -1), 
            0
        ) 
    );

    const volVectorField& cellCenters = mesh.C(); 
    const surfaceVectorField& faceCenters = mesh.Cf(); 

    vector circleCenter (0.5, 0.5, 0); 

    // Set internal face centered curvature field.
    forAll(faceCurvature, I)
    {
        faceCurvature[I] = circle_curvature(faceCenters[I], circleCenter);  
    }

    // TODO: set the boundary face centered curvature field

    
    // Set internal cell centered curvature field.
    forAll(cellCurvature, I)
    {
        cellCurvature[I] = circle_curvature(cellCenters[I], circleCenter);  
    }

    // TODO: set the boundary cell centered curvature field

    // Write the fields. 
    faceCurvature.write(); 
    cellCurvature.write();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
