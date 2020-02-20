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
#include "centroidIsoPointCalculator.H"
#include "linearLeastSquaresIsoPointCalculator.H"

// Test functions
// - Implicit representation
#include "sphereHypersurface.H"
#include "ellipsoidHypersurface.H"
// - Signed distance functions
#include "analyticalEllipsoid.H"
#include "analyticalPlane.H"
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

    OFstream centroidErrorFile(casePath + "/centroidPositioningErrors.csv"); 
    centroidErrorFile << "SURFACE,LINF_EDGE,L1_EDGE,L2_EDGE, LINF_CELL, L1_CELL, L2_CELL, CPU_TIME_SECONDS" << endl; 

    OFstream leastSquaresNoWeightingErrorFile(casePath + "/leastSquaresNoWeightingPositioningErrors.csv"); 
    leastSquaresNoWeightingErrorFile << "SURFACE,LINF_EDGE,L1_EDGE,L2_EDGE, LINF_CELL, L1_CELL, L2_CELL, CPU_TIME_SECONDS" << endl; 

    OFstream leastSquaresWeightedErrorFile(casePath + "/leastSquaresWeightedPositioningErrors.csv"); 
    leastSquaresWeightedErrorFile << "SURFACE,LINF_EDGE,L1_EDGE,L2_EDGE, LINF_CELL, L1_CELL, L2_CELL, CPU_TIME_SECONDS" << endl; 

    OFstream rbfErrorFile(casePath + "/rbfPositioningErrors.csv"); 
    rbfErrorFile << "RBF,STENCIL,SURFACE,LINF_CELL,L1_CELL,L2_CELL,POINT_CORR_CPU_TIME_SEC,FACTOR_CPU_TIME_SEC,SOL_CPU_TIME_SEC" << endl; 

    // Plane, Signed Distance
    surfaceTestFields fields(mesh, sigDistPlane);
    // Linear + Centroid reconstruction
    testIsoPoints<centroidIsoPointCalculator>(sigDistPlane, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sigDistPlane, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sigDistPlane, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF + Dual Contouring reconstruction
    {
        rbfReconstructLoop<rbfTuple>(sigDistPlane, fields, rbfErrorFile, casePath);
    }

    // Sphere, Signed Distance, Bulk 
    fields.setValues(sigDistSphere);
    // Linear + Centroid reconstruction
    testIsoPoints<centroidIsoPointCalculator>(sigDistSphere, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sigDistSphere, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sigDistSphere, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF + Dual Contouring reconstruction
    {
        rbfReconstructLoop<rbfTuple>(sigDistSphere, fields, rbfErrorFile, casePath);
    }

    // Sphere, Signed Distance, Boundary 
    fields.setValues(bSigDistSphere);
    // Linear + Centroid reconstruction
    testIsoPoints<centroidIsoPointCalculator>(bSigDistSphere, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bSigDistSphere, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bSigDistSphere, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF + Dual Contouring reconstruction
    {
        rbfReconstructLoop<rbfTuple>(bSigDistSphere, fields, rbfErrorFile, casePath); 
    }

    // Ellipsoid, Signed Distance, Bulk 
    fields.setValues(sigDistEllipsoid);
    // Linear + Centroid 
    testIsoPoints<centroidIsoPointCalculator>(sigDistEllipsoid, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sigDistEllipsoid, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(sigDistEllipsoid, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF
    {
        rbfReconstructLoop<rbfTuple>(sigDistEllipsoid, fields, rbfErrorFile, casePath);
    }

    // Ellipsoid, Signed Distance, Boundary 
    fields.setValues(bSigDistEllipsoid);
    // Linear + Centroid 
    testIsoPoints<centroidIsoPointCalculator>(bSigDistEllipsoid, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bSigDistEllipsoid, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bSigDistEllipsoid, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF
    {
        rbfReconstructLoop<rbfTuple>(bSigDistEllipsoid, fields, rbfErrorFile, casePath); 
    }

    // Sphere, Implicit, Bulk
    fields.setValues(implicitSphere);
    // Linear + Centroid
    testIsoPoints<centroidIsoPointCalculator>(implicitSphere, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(implicitSphere, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(implicitSphere, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF
    {
        rbfReconstructLoop<rbfTuple>(implicitSphere, fields, rbfErrorFile, casePath); 
    }

    // Sphere, Implicit, Boundary 
    fields.setValues(bImplicitSphere);
    // Linear + Centroid
    testIsoPoints<centroidIsoPointCalculator>(bImplicitSphere, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bImplicitSphere, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bImplicitSphere, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF
    {
        rbfReconstructLoop<rbfTuple>(bImplicitSphere, fields, rbfErrorFile, casePath); 
    }

    // Ellipsoid, Implicit, Bulk 
    fields.setValues(implicitEllipsoid);
    // Linear + Centroid
    testIsoPoints<centroidIsoPointCalculator>(implicitEllipsoid, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(implicitEllipsoid, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(implicitEllipsoid, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF
    {
        rbfReconstructLoop<rbfTuple>(implicitEllipsoid, fields, rbfErrorFile, casePath); 
    }

    // Ellipsoid, Implicit, Boundary 
    fields.setValues(bImplicitEllipsoid);
    // Linear + Centroid
    testIsoPoints<centroidIsoPointCalculator>(bImplicitEllipsoid, fields, centroidErrorFile, casePath); 
    // Linear + Linear Least Squares without weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bImplicitEllipsoid, fields, leastSquaresNoWeightingErrorFile, casePath, false);
    // Linear + Linear Least Squares with weighting
    testIsoPoints<linearLeastSquaresIsoPointCalculator>(bImplicitEllipsoid, fields, leastSquaresWeightedErrorFile, casePath, true);
    // RBF
    {
        rbfReconstructLoop<rbfTuple>(bImplicitEllipsoid, fields, rbfErrorFile, casePath); 
    }

    Info<< nl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
