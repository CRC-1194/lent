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

struct labelledPoints
{
    DynamicList<point> points; 
    DynamicList<label> labels; 

    void clear() { points.clear(); labels.clear(); };

    void append(const point& p, const label& l)
    {
        points.append(p); 
        labels.append(l); 
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
                  { return (n++ * (end - start)) / (result.size() - 1);}); 

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
    volScalarField& icells, // Visualization only.
    DynamicList<label> icellLabels
)
{
    icells = 0; 

    icellLabels.clear();

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
                icells[cellI] = 1;
                icellLabels.append(cellI); 
            }
        }
    }; 
}

void intersectedCells
(
    const volScalarField& dist, 
    const pointScalarField& pdist, 
    const DynamicList<label> iCellLabels,
    labelledPoints& iPointData 
)
{
    tmp<volVectorField> distGradTmp = fvc::grad(dist); 
    tmp<volScalarField> distGradDotTmp = distGradTmp() & distGradTmp(); 
    tmp<volVectorField> dirTmp = distGradTmp() / distGradDotTmp();  

    const volVectorField& dir = dirTmp(); 

    const fvMesh& mesh = dist.mesh(); 
    const volVectorField& C = mesh.C();  

    iPointData.clear(); 

    forAll(iCellLabels, iCell)
    {
        const label& iCellLabel = iCellLabels[iCell]; 
        const point& xc = C[iCellLabel]; 
        iPointData.append
        (
            xc - dist[iCellLabel] * dir[iCellLabel],
            iCellLabel
        );
    }
}

template<typename Points, typename String> 
void pointsToVtk
(
    Points const & points, 
    String const& fileName
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

std::string indexedString(const std::string& s, label index, label nPad=4)
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
    scalarField thetaSpace = linspace(0,pi,100); 
    scalarField phiSpace = linspace(0,2*pi,100); 

    // Test the plane reconstruction that is randomized in terms of the position 
    // in the unit box domain and with the parameterized orientation. 
    // The parameterization also contains normal orientations that are collinear 
    // with the coordinate axes.
    //point planePoint (dis(gen), dis(gen), dis(gen));
    point planePoint (0.5, 0.5, 0.5);

    DynamicList<label> icellList; 

    for (const auto& theta : thetaSpace)
    { 
        for (const auto& phi : phiSpace)
        {
            runTime++; 

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

            distance(dist, plane); 
            distance(pdist, plane);

            intersectedCells(dist, pdist, icells, icellList); 
            isoSurface isurf(dist, pdist, 0, false); 

            if (runTime.writeTime())
                isurf.write(indexedString("isosurf", runTime.timeIndex()) + ".stl"); 

            runTime.write(); 
        }
           
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
