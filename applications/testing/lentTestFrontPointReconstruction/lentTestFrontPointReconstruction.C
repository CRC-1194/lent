/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  OpenFOAM-plus 
    \\  /    A nd           | Copyright Tomislav Maric, TU Darmstadt 
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
#include "DynamicField.H"

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
    DynamicList<point> points_; 
    std::set<label> cellLabels_; 

    void clear() { points_.clear(); cellLabels_.clear(); };

    void clearPoints() { points_.clear(); }; 

    void clearLabels() { cellLabels_.clear(); };

    void append(const point& p, const label& l)
    {
        points_.append(p); 
        cellLabels_.insert(l); 
    }

    void append(const point& p)
    {
        points_.append(p); 
    }

    void append(const label& l)
    {
        cellLabels_.insert(l);
    }

    void append(const label& l, const point& p)
    {
        append(p,l); 
    }

    label nPoints() const
    {
        return points_.size(); 
    }

    label nCellLabels() const
    {
        return cellLabels_.size(); 
    }

    const DynamicList<point>& points() const
    {
        return points_; 
    }

    void resizePoints(const label& newsize)
    {
        points_.resize(newsize);
    }

    const std::set<label>& cellLabels() const
    {
        return cellLabels_;
    }

    void setPoint(const label& l, const point& p)
    {
        points_[l] = p; 
    }

    bool emptyLabels()
    {
        return cellLabels_.empty(); 
    }

    bool emptyPoints()
    {
        return points_.empty(); 
    }
};

struct surfErr
{
    scalar phi; 
    scalar theta; 
    scalar L1; 
    scalar L2; 
    scalar Linf;

    surfErr()
        : phi(0), theta(0), L1(0), L2(0), Linf(0)
    {};

    friend std::ostream& operator <<(std::ostream& os, const surfErr& errors); 
};

std::ostream& operator <<(std::ostream& os, const surfErr& errors)
{
    os << errors.phi << "," << errors.theta << "," 
        << errors.L1 << "," << errors.L2 << ","
        << errors.Linf << "\n";
    return os; 
} 

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

    forAll(pf.boundaryField(), patchI)
    {
        pointPatchField<scalar>& pfb = pf.boundaryFieldRef()[patchI]; 

        const pointPatch& patch = pfb.patch(); 
        const pointField& localPoints = patch.localPoints(); 
        const labelList& meshPoints = patch.meshPoints(); 
        forAll(pfb, pointI)
        {
            pf[meshPoints[pointI]] = surf.signedDistance(localPoints[pointI]);
        }
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
                iData.append(cellI); 
            }
        }
    }; 
}

void surfacePoints 
(
    const volScalarField& dist, 
    const pointScalarField& pdist, 
    surfaceData& iData 
)
{
    if (iData.emptyLabels())
        return; 

    // Clear only interface points and related data, as interface cell labels
    // are assumed to be computed and available.  
    
    // TODO: Check this: the loop below must loop over cellLabels and assign
    // each value, otherwise old point values may be left in the list. TM. 
    iData.resizePoints(iData.nCellLabels());

    tmp<volVectorField> distGradTmp = fvc::grad(dist); 
    tmp<volScalarField> distGradDotTmp = distGradTmp() & distGradTmp(); 
    tmp<volVectorField> dirTmp = distGradTmp() / (distGradDotTmp() + SMALL); 

    const volVectorField& dir = dirTmp(); 

    const fvMesh& mesh = dist.mesh(); 
    const volVectorField& C = mesh.C();  

    label labelI = 0; 
    for(const auto& cutCellI : iData.cellLabels())
    {
        iData.setPoint(labelI, C[cutCellI] - dist[cutCellI] * dir[cutCellI]);
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
        
void difference(const surfaceData& surf, const analyticalSurface& asurf, scalarField& diff) 
{
    const auto& surfPoints = surf.points(); 

    forAll(diff, pointI)
    {
        // The iso-surface value shold be 0 here: differencce
        // between the exact signed distance and the one given by the
        // directed derivative (0). TM.
        diff[pointI] = asurf.signedDistance(surfPoints[pointI]);
    }
}

int main(int argc, char **argv)
{
    argList::addOption
    (
        "Ntheta",
        "100",
        "Number of intervals for the spherical angle theta." 
    );

    argList::addOption
    (
        "Nphi",
        "100",
        "Number of intervals for the spherical angle phi." 
    );

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
    const label Ntheta = args.optionLookupOrDefault<label>("Ntheta", 100);
    const label Nphi = args.optionLookupOrDefault<label>("Nphi", 100);
    scalarField thetaSpace = linspace(0, pi, Ntheta); 
    scalarField phiSpace = linspace(0, 2*pi, Nphi); 

    // Test the plane reconstruction that is randomized in terms of the position 
    // in the unit box domain and with the parameterized orientation. 
    // The parameterization also contains normal orientations that are collinear 
    // with the coordinate axes.
    // TODO: Randomize the point position.
    //point planePoint (dis(gen), dis(gen), dis(gen));
    point planePoint (0.5, 0.5, 0.5);

    // Initialize surface data.
    surfaceData surfDat;

    // Open the error file for output and write the header line. 
    const label Nc = mesh.nCells(); 
    const scalar h = average(pow(mesh.deltaCoeffs(), -1)).value(); 
    OFstream errFile("dual-contouring-errors.dat"); 
    errFile << "Nc,h,phi,theta,Einf,E2" << nl;

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

            // Find cells cut by the surface.
            intersectedCells(dist, pdist, icells, surfDat); 
            
            // Compute interface points.
            surfacePoints(dist, pdist, surfDat);

            // Calculate the difference between the reconstructed surface
            // and the analytical surface. 
            scalarField diff(surfDat.nPoints(), 0); 
            difference(surfDat, plane, diff); 

            // Calculate and write the error values. 
            const scalar Linf = max(mag(diff));  
            const scalar L2 = Foam::sqrt(sum(sqr(diff)));
            errFile << Nc << "," << h << "," << phi 
                << "," << theta << ","  << Linf << "," << L2 << nl; 

            // Write selected test data.
            if (runTime.writeTime() || runTime.timeIndex() == 0)
            {
                isurf.write(indexedString("isosurf", runTime.timeIndex()) + ".stl"); 
                icells.write();
                dist.write(); 
                pdist.write(); 
                writePointsToVtk(
                    surfDat.points(),
                    indexedString("points", runTime.timeIndex()) + ".vtk"
                );

            }
            runTime++; 
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
