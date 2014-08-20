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

scalar circle_curvature(const point& p, const point& circleCenter)
{
    vector distance = circleCenter - p;

    // Since the circle is used in 2D test cases, project the vector
    // on the x-y plane
    distance.z() = 0;

    double radius = Foam::sqrt(distance & distance);
    double curvature = 0;

    // To avoid divide y zero in case the circle center and the point
    // coincide, check distance and limit the curvature value
    //
    // TODO: think of reasonable limits for the minimal radius and
    // maximum curvature. Maybe the value for minimal radius can be based
    // on the edge length of the cells.
    if (radius > 1.0e-10) curvature = 1 / radius;
    else                  curvature = 1.0e10;

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
        "x",
        "x-coordinate of the shape center"
    );   

    argList::addOption
    (
        "y",
        "y-coordinate of the shape center"
    );

    argList::addOption
    (
        "z",
        "z-coordinate of the shape center"
    );   

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (!args.optionFound("shape"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-shape' to select the shape for which the "
            << "exact is to be calculated." << endl << exit(FatalError);
    }

    if (!args.optionFound("x"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-x' to specify the x-coordinate of the"
            << " shape center." << endl << exit(FatalError);
    }

    if (!args.optionFound("y"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-y' to specify the y-coordinate of the"
            << " shape center." << endl << exit(FatalError);
    }

    if (!args.optionFound("z"))
    {
        FatalErrorIn
        (
            "main()"
        )   << "Please use option '-z' to specify the z-coordinate of the"
            << " shape center." << endl << exit(FatalError);
    }

    const word shapeName = args.optionRead<word>("shape");
    const double x_center = args.optionRead<double>("x");
    const double y_center = args.optionRead<double>("y");
    const double z_center = args.optionRead<double>("z");

    vector circleCenter (x_center, y_center, z_center);

    // Get the time directories from the simulation folder using time selector
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        // Set the time to the current time directory and set the time index
        // to timeI
        runTime.setTime(timeDirs[timeI], timeI);

        // Initialize the field objects for curvature
        surfaceScalarField faceCurvatureExact  
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

        // Debugging/Testing: write coordinates of face centers to a field
        surfaceVectorField faceCenterCoordinates
        (
            IOobject
            (
                "faceCenterCoordinates",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector(
                "zero", 
                dimensionSet(1,0,0,0,0,0,0),
                Foam::vector(0,0,0)
            ) 
        );

        volScalarField cellCurvatureExact  
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

        const volVectorField& cellCenters = mesh.C(); 
        const surfaceVectorField& faceCenters = mesh.Cf(); 

        // Set internal face centered curvature field.
        forAll(faceCurvatureExact, I)
        {
            faceCurvatureExact[I] = circle_curvature(faceCenters[I],
                                                     circleCenter);
            // Debugging
            faceCenterCoordinates[I] = faceCenters[I];
        }

        // TODO: set the boundary face centered curvature field
        // In this loop something is rotten. Probably the inner loop takes
        // the wrong face center coordinates
        forAll(faceCurvatureExact.boundaryField(), I)
        {
            scalarField& bScalarField = faceCurvatureExact.boundaryField()[I];

            // Debugging
            vectorField& bVectorField = faceCenterCoordinates.boundaryField()[I];

            forAll(bScalarField, K)
            {  
                bScalarField[K] = circle_curvature(faceCenters[K],
                                                   circleCenter);
                // Debugging
                bVectorField[K] = faceCenters[K];
            }

            /*
            // Another approach, yields a compile-time error when executing
            // wmake
            
            const fvMesh& mesh = bScalarField.mesh();
            const surfaceVectorField& boundaryFaceCenters = mesh.Cf();
            
            forAll(bScalarField, K)
            {
                bScalarField[K] = circle_curvature(boundaryFaceCenters[K],
                                                   circleCenter);
                // Debugging
                bVectorField[K] = boundaryFaceCenters[K];
            }
            */
        }
    
        // Set internal cell centered curvature field.
        forAll(cellCurvatureExact, I)
        {
            cellCurvatureExact[I] = circle_curvature(cellCenters[I],
                                                     circleCenter);
        }

        // Set the boundary cell centered curvature field
        forAll(cellCurvatureExact.boundaryField(), I)
        {
            scalarField& bScalarField = cellCurvatureExact.boundaryField()[I];

            forAll(bScalarField, I)
            {
                bScalarField[I] = circle_curvature(cellCenters[I],
                                                   circleCenter);
            }
        }

        // Write the fields. 
        faceCurvatureExact.write(); 
        cellCurvatureExact.write();

        // Debugging
        faceCenterCoordinates.write();
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
