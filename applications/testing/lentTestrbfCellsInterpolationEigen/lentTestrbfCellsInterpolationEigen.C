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

#include "lentTestrbfCellsInterpolationEigenTemplates.H"

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
            << "where N is the size of the RBF system.\n";

    #include "createFields.H"

    const auto& meshPoints = mesh.points(); 
    const auto& cellCenters = mesh.C(); 

    // Error file.
    OFstream errorFile ("lentTestrbfCellsInterpolationEigen.dat"); 
    errorFile << "RBF,SURFACE,LINF_BCC,LINF_BCC-FVM" << endl; 

    using rbfKernelType = std::tuple_element_t<0, rbfTuple>;
    
    Info << "RBF Kernel = " << rbfKernelType::name() << endl;

    // Initialize and factorize the RBF interpolation linear equation systems. 
    Info << "Factorizing the interpolation matrices... "; 
    rbfCellsInterpolationEigen<rbfKernelType> cellRbfsBcc(mesh, stencilType::BCC); 
    rbfCellsInterpolationEigen<rbfKernelType> cellRbfsBccFvm(mesh, stencilType::BCC_FVM); 
    Info << "done." << endl; 

    for (decltype(radii.size()) testI = 0; testI < radii.size(); ++testI)
    {
        runTime.setTime(testI, testI); 
        Info<< "Test = " << runTime.timeName() << nl << endl;

        // Reset the errors.
        LinfEllipsoidBcc == linfInit;  
        LinfEllipsoidBccFvm == linfInit;
        LinfSphereBcc == linfInit;  
        LinfSphereBccFvm == linfInit;

        // BEGIN ELLIPSOID TEST
        // - Set signed distance fields.
        ellipsoidHypersurface ellipsoid (aAxes[testI], bAxes[testI], cAxes[testI]); 
        surfaceSetField(veField, cellCenters, ellipsoid);
        surfaceSetField(peField, meshPoints, ellipsoid);
        
        // - Interpolate new field values. 
        cellRbfsBcc.solve(veField, peField); 
        evaluateLinfErrors(LinfEllipsoidBcc, randPointsInCells, 
                           ellipsoid, cellRbfsBcc);

        cellRbfsBccFvm.solve(veField, peField); 
        evaluateLinfErrors(LinfEllipsoidBccFvm, randPointsInCells, 
                           ellipsoid, cellRbfsBccFvm);

        // Report ellipsoid errors.
        auto LinfBcc = max(LinfEllipsoidBcc).value(); 
        auto LinfBccFvm = max(LinfEllipsoidBccFvm).value(); 
        errorFile << rbfKernelType::name() << "," << "ELLIPSOID," 
            << LinfBcc << "," << LinfBccFvm << endl;
        Info << "Linf ellipsoid, BCC stencil = " 
            << LinfBcc << endl;
        Info << "Linf ellipsoid, BCC_FVM stencil = " 
            << LinfBccFvm << endl;

        // END ELLIPSOID TEST

        // BEGIN SPHERE TEST
        // - Set signed distance fields.
        sphereHypersurface sphere(point(0.,0.,0.), radii[testI]); 
        surfaceSetField(psField, meshPoints, sphere);
        surfaceSetField(vsField, cellCenters, sphere);

        // - Interpolate new field values. 
        cellRbfsBcc.solve(vsField, psField); 
        evaluateLinfErrors(LinfSphereBcc, randPointsInCells, 
                           sphere, cellRbfsBcc);

        cellRbfsBccFvm.solve(vsField, psField); 
        evaluateLinfErrors(LinfSphereBccFvm, randPointsInCells, 
                           sphere, cellRbfsBccFvm);
        
        // Report sphere errors.
        LinfBcc = max(LinfSphereBcc).value(); 
        LinfBccFvm = max(LinfSphereBccFvm).value(); 
        errorFile << rbfKernelType::name() << "," << "SPHERE," 
            << LinfBcc << "," << LinfBccFvm << endl;
        Info << "Linf sphere, BCC stencil = " << LinfBcc << endl;
        Info << "Linf sphere, BCC_FVM stencil = " << LinfBccFvm << endl;
        // END SPHERE TEST
        
        runTime.writeNow(); 
        runTime.printExecutionTime(Info);
    }

    Info<< nl;
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
