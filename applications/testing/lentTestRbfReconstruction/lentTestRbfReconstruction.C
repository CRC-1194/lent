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
#include "rbfIsoPointCalculator.H"
#include "isoPointCalculator.H"

// Test functions
// - Implicit representation
#include "sphereHypersurface.H"
#include "ellipsoidHypersurface.H"
// - Signed distance functions
#include "analyticalEllipsoid.H"
#include "analyticalSphere.H"

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

    #include "testSurfaces.H"

    std::string casePath = args.rootPath() + "/" + args.globalCaseName();

    // TODO: Rename the isoPointCalculator to centroidIsoPointCalculator. TM.
    // Test linear point reconstruction: 
    OFstream centroidErrorFile(casePath + "/centroidPositioningErrors.csv"); 
    centroidErrorFile << "SURFACE,LINF_EDGE,L1_EDGE,L2_EDGE, LINF_CELL, L1_CELL, L2_CELL, CPU_TIME_SECONDS" << endl; 

    // Signed distance sphere test: bulk and boundary  
    testIsoPoints(sigDistSphere, sigDistSphereFields, centroidErrorFile, casePath); 
    testIsoPoints(bSigDistSphere, bSigDistSphereFields, centroidErrorFile, casePath); 

    // Signed distance ellipsoid test: bulk and boundary 
    testIsoPoints(sigDistEllipsoid, sigDistEllipsoidFields, centroidErrorFile, casePath); 
    testIsoPoints(bSigDistEllipsoid, bSigDistEllipsoidFields, centroidErrorFile, casePath); 

    // Implicit sphere test: bulk and boundary 
    testIsoPoints(implicitSphere, implicitSphereFields, centroidErrorFile, casePath); 
    testIsoPoints(bImplicitSphere, bImplicitSphereFields, centroidErrorFile, casePath); 

    // Implicit ellipsoid test: bulk and boundary 
    testIsoPoints(implicitEllipsoid, implicitEllipsoidFields, centroidErrorFile, casePath); 
    testIsoPoints(bImplicitEllipsoid, bImplicitEllipsoidFields, centroidErrorFile, casePath); 
    
    OFstream rbfErrorFile(casePath + "/rbfPositioningErrors.csv"); 
    rbfErrorFile << "RBF,STENCIL,SURFACE,LINF_CELL,L1_CELL,L2_CELL,POINT_CORR_CPU_TIME_SEC,FACTOR_CPU_TIME_SEC,SOL_CPU_TIME_SEC" << endl; 

    // Test RBF reconstruction: loop over all RBF kernels at compile time.
    
    // Signed distance sphere test: bulk and boundary  
    rbfReconstructLoop<rbfTuple>(sigDistSphere, sigDistSphereFields, rbfErrorFile, casePath); 
    rbfReconstructLoop<rbfTuple>(bSigDistSphere, bSigDistSphereFields, rbfErrorFile, casePath); 

    // Signed distance ellipsoid test: bulk and boundary 
    rbfReconstructLoop<rbfTuple>(sigDistEllipsoid, sigDistEllipsoidFields, rbfErrorFile, casePath); 
    rbfReconstructLoop<rbfTuple>(bSigDistEllipsoid, bSigDistEllipsoidFields, rbfErrorFile, casePath); 

    // Implicit sphere test: bulk and boundary 
    rbfReconstructLoop<rbfTuple>(implicitSphere, implicitSphereFields, rbfErrorFile, casePath); 
    rbfReconstructLoop<rbfTuple>(bImplicitSphere, bImplicitSphereFields, rbfErrorFile, casePath); 

    // Implicit ellipsoid test: bulk and boundary 
    rbfReconstructLoop<rbfTuple>(implicitEllipsoid, implicitEllipsoidFields, rbfErrorFile, casePath); 
    rbfReconstructLoop<rbfTuple>(bImplicitEllipsoid, bImplicitEllipsoidFields, rbfErrorFile, casePath); 

    Info<< nl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
