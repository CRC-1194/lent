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
#include "centroidIsoPointCalculator.H"
#include "linearLeastSquaresIsoPointCalculator.H"

// Time measurement.
#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
using namespace RBF;
using namespace FrontTracking;

#include "lentTestRbfReconstructionHelpers.H"

int main(int argc, char **argv)
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    if(mesh.nCells() > 64 * 64 * 64) 
        WarningIn(__PRETTY_FUNCTION__) 
            << "Request for memory may be too large. " << endl
            << "Memory usage of the RBF factorization is NxN x nCells x 64 bits, " << endl
            << "where N is the number of points in the RBF stencil (hex cells: 9 for BCC, 15 for BCC_FVM).\n";

    // Initialize the surface test fields.
    ellipsoidHypersurface ellipsoid (1/3., 1/2., 2/3.);
    surfaceTestFields ellipsoidFields(mesh, ellipsoid, "Ellipsoid");
    sphereHypersurface sphere(point(0., 0., 0.), 0.5);
    surfaceTestFields sphereFields(mesh, sphere, "Sphere");

    std::string casePath = args.rootPath() + "/" + args.globalCaseName();

    // Test centroid point reconstruction: 
    OFstream centroidErrorFile(casePath + "/centroidPositioningErrors.csv"); 
    centroidErrorFile << "SURFACE,LINF_EDGE,L1_EDGE,L2_EDGE, LINF_CELL, L1_CELL, L2_CELL, CPU_TIME_SECONDS" << endl; 
    testIsoPoints<centroidIsoPointCalculator>(sphere, sphereFields, centroidErrorFile, casePath, false); 
    testIsoPoints<centroidIsoPointCalculator>(ellipsoid, ellipsoidFields, centroidErrorFile, casePath, false); 

    // Test linear least squares point reconstruction without weighting
    OFstream leastSquaresErrorFile(casePath + "/leastSquaresNoWeightingPositioningErrors.csv"); 
    leastSquaresErrorFile << "SURFACE,LINF_EDGE,L1_EDGE,L2_EDGE, LINF_CELL, L1_CELL, L2_CELL, CPU_TIME_SECONDS" << endl; 
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sphere, sphereFields, leastSquaresErrorFile, casePath, false); 
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(ellipsoid, ellipsoidFields, leastSquaresErrorFile, casePath, false); 
    
    // Test linear least squares point reconstruction 
    OFstream leastSquaresWeightedErrorFile(casePath + "/leastSquaresWeightedPositioningErrors.csv"); 
    leastSquaresWeightedErrorFile << "SURFACE,LINF_EDGE,L1_EDGE,L2_EDGE, LINF_CELL, L1_CELL, L2_CELL, CPU_TIME_SECONDS" << endl; 
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sphere, sphereFields, leastSquaresWeightedErrorFile, casePath, true); 
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(ellipsoid, ellipsoidFields, leastSquaresWeightedErrorFile, casePath, true); 


    OFstream rbfErrorFile(casePath + "/rbfPositioningErrors.csv"); 
    rbfErrorFile << "RBF,STENCIL,SURFACE,LINF_CELL,L1_CELL,L2_CELL,POINT_CORR_CPU_TIME_SEC,FACTOR_CPU_TIME_SEC,SOL_CPU_TIME_SEC" << endl; 

    // Test RBF reconstruction: loop over all RBF kernels at compile time.
    rbfReconstructLoop<rbfTuple>(sphere, sphereFields, rbfErrorFile, casePath); 
    rbfReconstructLoop<rbfTuple>(ellipsoid, ellipsoidFields, rbfErrorFile, casePath); 

    ellipsoidFields.writeValueFields();
    sphereFields.writeValueFields();

    Info<< nl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
