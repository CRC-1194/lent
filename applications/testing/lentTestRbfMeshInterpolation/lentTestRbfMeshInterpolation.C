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
#include "analyticalEllipsoid.H"
#include "sphereHypersurface.H"
#include "ellipsoidHypersurface.H"
#include "rbfIsoPointCalculator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace RBF;
using namespace FrontTracking;

#include "lentTestRbfMeshInterpolationTemplates.H"

int main(int argc, char **argv)
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if(mesh.nCells() > 64 * 64 * 64) 
        WarningIn(__PRETTY_FUNCTION__) 
            << "Request for memory may be too large. " << endl
            << "The RBF system factorization is memory intensive," 
            << "the program may run out of memory." << endl
            << "Memory usage of the RBF factorization is NxN x nCells x 64 bits, " << endl
            << "where N is the number of points in the RBF stencil.\n";

    // Initialize the surface test fields.
    surfaceTestFields ellFields(mesh, "Ellipsoid"); 
    surfaceTestFields sphFields(mesh, "Sphere"); 

    // Populate mesh cells with nPointsPerCell random points.
    const label nPointsPerCell = 10; 
    Info << "Generating " << nPointsPerCell << " random sampling points per cell ..."; 
    using pointVectorVector = std::vector<std::vector<point>>;
    pointVectorVector randPointsInCells
    (
        mesh.nCells(), 
        std::vector<point>(nPointsPerCell, point(VGREAT, VGREAT, VGREAT))
    );
    genRandomPointsInCells(randPointsInCells, mesh, nPointsPerCell);
    Info << "Done." << endl;

    // Error file.
    OFstream errorFile ("lentTestrbfCellsInterpolationEigen.dat"); 
    errorFile << "RBF,SURFACE,LINF_BCC,LINF_BCC-FVM" << endl; 

    // Tests the RBF interpolation with an ellipsoid and sphere hypersurface
    // whose radii increase with respect to the mesh size. 
    testRbfEllipsoidSphere<rbfTuple>
    (
        mesh, 
        runTime,  
        randPointsInCells,
        ellFields,
        sphFields,
        errorFile
    );

    Info<< nl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
