/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | Copyright (C) 2018 Tomislav Maric, TU Darmstadt
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
    Tomislav Maric maric@mma.tu-darmstadt.de

Description
    Test application for the RBF interpolation on an FV mesh.

    Tests
        - rbfCellsInterpolationEigen: RBF interpolation over cell stencils

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"

#include "rbFunctions.H"
#include "rbfCellsInterpolationEigen.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace RBF;

int main(int argc, char **argv)
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Cell centered (volume) field.
    volScalarField cellField 
    (
        IOobject
        (
            "cellField",
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ), 
        // FIXME: Replace this with a function.
        Foam::pow(mesh.C().component(0),2) + 
        Foam::pow(mesh.C().component(1),2) + 
        Foam::pow(mesh.C().component(2),2) 
    );

    // Point-centered (cell-corner) field.
    pointMesh pointMesh(mesh); 
    pointScalarField cornerField 
    (
        IOobject
        (
            "cornerField",
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ), 
        pointMesh, 
        dimensionedScalar("cornerField", Foam::pow(dimLength,2), 0)
    );

    // Write the fields.
    cornerField.write();
    cellField.write(); 

    // FIXME: Abstract this into some kind of a setter function for the point field. 
    const auto& meshPoints = mesh.points(); 
    forAll(cornerField, pointI)
    {
        cornerField[pointI] = Foam::pow(meshPoints[pointI][0], 2) +
            Foam::pow(meshPoints[pointI][1], 2) +
            Foam::pow(meshPoints[pointI][1], 2);
    }

    // For all kernels...
    const std::string rbfName ("MULTIQUADRIC");
    auto rbfKernel = [rbfName](double r, double rs) { return rbfKernels.at(rbfName)(r, rs, 1.0); };  

    // Construct 
    // - Assemble the point stencils from the mesh and set the RBF interpolation kernels.
    // - Factorize the lineary systems based on the set kernel and calculated stencils. 
    rbfCellsInterpolationEigen cellRbfs(mesh, rbfKernel, stencilType::BCC_FVM); 
    //cellRbfs.interpolate(cellField, cornerField); 


    // TESTS: 
    // For each cell, sample the RBF interpolant randomly, evaluate Einf and Ermse
    //      - Store Einf and Ermse fields for visualization.
    //      - Compute the E_inf for the mesh. 
    //      - Write error data into a file: Ncells, E_inf 

    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
