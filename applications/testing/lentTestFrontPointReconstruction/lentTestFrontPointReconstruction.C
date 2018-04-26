/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  OpenFOAM-plus 
    \\  /    A nd           | Copyright held by original author
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
    Test different reconstruction algorithms for computing front points from 
    the signed distance.

\*---------------------------------------------------------------------------*/

// Algorithm
#include "fvCFD.H"

#include "pointFields.H"
#include "isoSurface.H"
#include <set>
 // - For storing interface points.
#include "DynamicList.H"

// Testing 
#include <random> 
#include <algorithm>
#include "mathematicalConstants.H"
#include "analyticalPlane.H"
#include "analyticalSphere.H"
#include <sstream>
#include <iomanip>
 //- VTK output 
#include "foamVtkLegacyAsciiFormatter.H"
#include "foamVtkOutput.H"
#include "pointList.H"
#include "point.H"
#include "List.H"
#include "pointList.H"
#include <fstream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO: Make this a class.
struct surfaceData 
{
    DynamicList<point> points; 
    // TODO: Use a set for the cell labels to avoid duplicate cells that are 
    //       generated when the interface passes exaclty through a mesh vertex.
    //DynamicList<label> cellLabels; 
    std::set<label> cellLabels; 

    void clear() { points.clear(); cellLabels.clear(); };

    void clearPoints() { points.clear(); }; 

    void clearLabels() { cellLabels.clear(); };

    void append(const point& p, const label& l)
    {
        points.append(p); 
        cellLabels.insert(l); 
    }

    void append(const label& l, const point& p)
    {
        append(p,l); 
    }
};

scalarField linspace(scalar start, scalar end, label size)
{
    scalarField result(size, 0); 

    std::generate(result.begin(), result.end(), 
                  [n = 0, start, end, result] () mutable 
                  { return (n++ * (end - start)) / 
                           (result.size() - 1 > 0 ? result.size() - 1 : 1);}); 

    return result; 
}

using namespace FrontTracking; 

void distance(volScalarField& vf, const analyticalSurface& surf)
{
    const auto& C = vf.mesh().C(); 

    forAll(vf, cellI)
    {
        vf[cellI] = surf.signedDistance(C[cellI]); 
    }

    const auto& Cf = vf.mesh().Cf(); 

    forAll(vf.boundaryField(), patchI)
    {
        fvPatchScalarField& vfb = vf.boundaryFieldRef()[patchI]; 
        const auto& patch = vfb.patch(); 

        forAll(vfb, faceI)
        {
            vfb[faceI] = surf.signedDistance(Cf[patch.start() + faceI]); 
        }
    }

}

void distance(pointScalarField& pf, const analyticalSurface& surf)
{
    const auto& points = pf.mesh()().points(); 

    forAll(pf, pointI)
    {
        pf[pointI] = surf.signedDistance(points[pointI]); 
    }
}

void intersectedCells
(
    const volScalarField& dist, 
    const pointScalarField& pdist, 
    volScalarField& iCells, // FIXME: Debugging, remove.
    surfaceData& iData 
)
{
    // Reset the intersected cells marker.
    iCells = 0; 
    // Clear everythng when new intersected cells are computed.
    iData.clear();

    const auto& mesh = dist.mesh(); 
    const auto& cellPointsLabelList = mesh.cellPoints(); 


    forAll(dist, cellI)
    {
        const auto& cellPointLabels = cellPointsLabelList[cellI];  

        const scalar& d0 = pdist[cellPointLabels[0]]; 

        forAll(cellPointLabels, pointI)
        {
            if (d0 * pdist[cellPointLabels[pointI]] <= SMALL)
            {
                iCells[cellI] = 1;
                iData.cellLabels.insert(cellI); 
            }
        }
    }; 
}

void interfacePoints 
(
    const volScalarField& dist, 
    const pointScalarField& pdist, 
    surfaceData& iData 
)
{
    if (iData.cellLabels.empty())
        return; 

    // Clear only interface points and related data, as interface cell labels
    // are assumed to be computed and available.  
    
    // TODO: Check this: the loop below must loop over cellLabels and assign
    // each value, otherwise old point values may be left in the list. TM. 
    iData.points.resize(iData.cellLabels.size());

    tmp<volVectorField> distGradTmp = fvc::grad(dist); 
    tmp<volScalarField> distGradDotTmp = distGradTmp() & distGradTmp(); 
    tmp<volVectorField> dirTmp = distGradTmp() / (distGradDotTmp() + SMALL); 

    const volVectorField& dir = dirTmp(); 

    const fvMesh& mesh = dist.mesh(); 
    const volVectorField& C = mesh.C();  

    label labelI = 0; 
    for(const auto& cutCellI : iData.cellLabels)
    {
        iData.points[labelI] = C[cutCellI] - dist[cutCellI] * dir[cutCellI];
        ++labelI;
    }
}

void writePointsToVtk
(
    UList<point> const & points, 
    std::string const& fileName
)
{
    // Open file for output in overwrite mode.
    std::ofstream vtkFile(fileName);
 
    // Initialize the VTK formatter, legacy in this case.
    vtk::legacyAsciiFormatter legacyFormat(vtkFile, 15);
 
    // Write the file header based on the chosen format (legacy VTK).
    vtk::legacy::fileHeader(legacyFormat, "points", vtk::fileTag::POLY_DATA);
 
    // Write the points beginning line for legacy VTK. 
    vtk::legacy::beginPoints(vtkFile, points.size());
 
    // Write the list of points using the legacy formatter.
    vtk::writeList(legacyFormat, points); 
}

std::string indexedString(const std::string& s, label index, label nPad=8)
{
    std::stringstream ss; 
    ss << s << "-" << std::setw(nPad) << std::setfill('0') << index;  
    return ss.str(); 
}

int main(int argc, char **argv)
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField dist 
    (
        IOobject
        (
            "dist", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ), 
        mesh
    );

    volScalarField icells 
    (
        IOobject
        (
            "icells", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ), 
        mesh
    );

    pointMesh pmesh(mesh); 
    
    pointScalarField pdist 
    (
        IOobject
        (
            "pdist", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ), 
        pmesh 
    );

    // Randomize plane position.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    // Parameterize plane orientation using spherical angles.
    const scalar pi = constant::mathematical::pi; 
    // TODO: Add options for nTheta and nPhi.
    scalarField thetaSpace = linspace(0,pi,10); 
    scalarField phiSpace = linspace(0,2*pi,10); 

    // Test the plane reconstruction that is randomized in terms of the position 
    // in the unit box domain and with the parameterized orientation. 
    // The parameterization also contains normal orientations that are collinear 
    // with the coordinate axes.
    //point planePoint (dis(gen), dis(gen), dis(gen));
    point planePoint (0.5, 0.5, 0.5);

    surfaceData surfDat;

    for (const auto& phi: thetaSpace)
    { 
        for (const auto& theta: phiSpace)
        {

            Info << "Iteration = " << runTime.timeIndex() << nl;

            analyticalPlane plane
            (
                planePoint, 
                vector
                (
                    Foam::sin(theta)*Foam::cos(phi), 
                    Foam::sin(theta)*Foam::sin(phi), 
                    Foam::cos(theta)
                )
            );

            // Set the cell and point distance fields using 
            // the analytical plane. 
            distance(dist, plane); 
            distance(pdist, plane);

            // Test the OpenFOAM iso-surface algorithm. 
            isoSurface isurf(dist, pdist, 0, false); 

            // Reset surface data.
            surfDat.clear(); 

            // Find cells cut by the surface.
            intersectedCells(dist, pdist, icells, surfDat); 
            
            // Compute interface points.
            interfacePoints(dist, pdist, surfDat);

            // Write selected test data.
            if (runTime.writeTime() || runTime.timeIndex() == 0)
            {
                isurf.write(indexedString("isosurf", runTime.timeIndex()) + ".stl"); 
                icells.write();
                dist.write(); 
                pdist.write(); 
                writePointsToVtk(
                    surfDat.points,
                    indexedString("points", runTime.timeIndex()) + ".vtk"
                );
                // TODO: Remove, testing.
                volVectorField dir ("dir", fvc::grad(dist) / ((fvc::grad(dist) & fvc::grad(dist)) + SMALL)); 
                dir.write(); 
            }

            runTime++; 
        }
           
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
