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
#include "timeSelector.H"

scalar circle_curvature(const point& p, const point& center)
{
    scalar curvature = 0;

    scalar radius = mag(center - p); 

    if (radius > SMALL) 
    {
        curvature = 1 / radius; 
    }
    else
    {
        curvature = GREAT; 
    }

    return curvature;
}

scalar sphere_curvature(const point& p, const point& sphereCenter)
{
    vector distance = sphereCenter - p;

    double radius = Foam::sqrt(distance & distance);
    double curvature = 0;

    // To avoid divide y zero in case the sphere center and the point
    // coincide, check distance and limit the curvature value
    //
    // TODO: think of reasonable limits for the minimal radius and
    // maximum curvature. Maybe the value for minimal radius can be based
    // on the edge length of the cells.
    if (radius > 1.0e-10) curvature = 2 / radius;
    else                  curvature = 2.0e10;

    return curvature;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    // Add necessary options and check if they have been set by the user
    argList::addOption
    (
        "shape",
        "Shape for which analytic curvature is to be set"
    );

    argList::addOption
    (
        "center",
        "Shape center position."
    );   


    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.optionFound("shape"))
    {
        FatalErrorIn("main()")   
            << "Please use option '-shape' to select the shape for which the "
            << "exact is to be calculated." << endl << exit(FatalError);
    }

    if (!args.optionFound("center"))
    {
        FatalErrorIn("main()")  
            << "Please use option '-center' to specify the shape center." 
            << endl << exit(FatalError);
    }

    const word shapeName = args.optionRead<word>("shape");

    vector center = args.optionRead<vector>("center"); 
    
    // Get the time directories from the simulation folder using time selector
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    const volVectorField& C = mesh.C(); 
    const surfaceVectorField& Cf = mesh.Cf(); 

    // Initialize the face centered curvature field.
    surfaceScalarField faceCurvature  
    (
        IOobject
        (
            "faceCurvatureExact",
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

    // Initialize the cell centered curvature field.
    volScalarField cellCurvature  
    (
        IOobject
        (
            "cellCurvatureExact",
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

    forAll(timeDirs, timeI)
    {
        // Set the time to the current time directory and set the time index to
        // timeI.
        runTime.setTime(timeDirs[timeI], timeI);

        // Set internal face centered curvature field.
        forAll(faceCurvature, I)
        {
            faceCurvature[I] = circle_curvature(Cf[I], center);
        }

        forAll(faceCurvature.boundaryField(), I)
        {
            fvsPatchScalarField& faceCurvatureBoundary = faceCurvature.boundaryField()[I];

            const fvsPatchVectorField& CfBoundary = Cf.boundaryField()[I];

            forAll(faceCurvatureBoundary, J)
            {  
                // FIXME: this was a bug, you have to address the CfBoundary
                //bScalarField[K] = circle_curvature(Cf[K], center);
                faceCurvatureBoundary[J] = circle_curvature(CfBoundary[J], center);
            }
        }
    
        // Set internal cell centered curvature field.
        forAll(cellCurvature, I)
        {
            cellCurvature[I] = circle_curvature(C[I], center);
        }

        // Set the boundary cell centered curvature field
        forAll(cellCurvature.boundaryField(), I)
        {
            fvPatchScalarField& cellCurvatureBoundary = cellCurvature.boundaryField()[I];
            const fvsPatchVectorField& CfBoundary = Cf.boundaryField()[I];

            forAll(cellCurvatureBoundary, J)
            {
                cellCurvatureBoundary[J] = circle_curvature(CfBoundary[J], center);
            }
        }

        // Write the fields. 
        faceCurvature.write(); 
        cellCurvature.write();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
