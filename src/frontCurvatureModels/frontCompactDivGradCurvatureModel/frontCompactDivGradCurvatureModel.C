/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Version:  2.2.x                               
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

Class
    Foam::frontCompactDivGradCurvatureModel

SourceFiles
    frontCurvatureModel.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description
    Compact version of the divGrad curvature Model. Uses the front-mesh
    communication to transfer the curvature values of interface cells
    to all cells in the narrow band

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#include "frontCompactDivGradCurvatureModel.H"

#include <algorithm>
#include <iterator>

#include "addToRunTimeSelectionTable.H"
#include "triSurfaceFields.H"

#include "lentCommunication.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontCompactDivGradCurvatureModel, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontCompactDivGradCurvatureModel, Dictionary);


// * * * * * * * * * * * * *  Private member functions * * * * * * * * * * * //
void frontCompactDivGradCurvatureModel::applySphereCorrection(
    const fvMesh& mesh,
    const triSurfaceFront& front,
    const volScalarField& signedDistance,
    volScalarField& curvature
) const
{
    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(front, mesh)
                                ); 
    // Only used as a list of all cells intersected by the interface
    const auto& interfaceCells = communication.interfaceCellToTriangles();

    for (const auto& p : interfaceCells)
    {
        const auto& temp = curvature[p.first];
        curvature[p.first] = sign(temp) * 2.0 /
                               (2.0/mag(temp) - signedDistance[p.first]);
    }
}

void frontCompactDivGradCurvatureModel::applyInterpolationCorrection(
    const fvMesh& mesh,
    const triSurfaceFront& front,
    const volScalarField& signedDistance,
    volScalarField& curvature
) const
{
    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(front, mesh)
                                ); 
    // Only used as a list of all cells intersected by the interface
    const auto& interfaceCells = communication.interfaceCellToTriangles();
    const auto& cellsTriangleNearest = communication.cellsTriangleNearest();
    const auto& C = mesh.C();

    for (const auto& p : interfaceCells)
    {
        const auto& cellLabel = p.first;

        // Idea: find neighbour cell for which the connection of the
        // two cell centres deviates the least from the interface normal
        //
        // Note: it is not sufficient to look only at the cells connected by a
        // face
        scalar maxCosine = 0.0;
        label interpolationCellLabel = -1;
        auto normal = cellsTriangleNearest[cellLabel].hitPoint() - C[cellLabel];
        normal /= (mag(normal) + SMALL);

        // No need to correct if cell centre is located very close to the
        // interface
        if (mag(signedDistance[cellLabel]) > SMALL)
        {
            auto extendedNeighbourCells = findNeighbourCells(mesh, cellLabel);

            for (const auto& neighbourCellLabel : extendedNeighbourCells)
            {
                if (signedDistance[cellLabel]*signedDistance[neighbourCellLabel] > 0.0)
                {
                    continue;
                }

                auto cc = C[cellLabel] - C[neighbourCellLabel];
                cc /= mag(cc);

                auto cosine = mag(cc & normal);

                if (cosine > maxCosine)
                {
                    maxCosine = cosine;
                    interpolationCellLabel = neighbourCellLabel;
                }
            }

            curvature[cellLabel] = 
                (
                 curvature[cellLabel]*mag(signedDistance[interpolationCellLabel])
                 +
                 curvature[interpolationCellLabel]*mag(signedDistance[cellLabel])
                )
                /
                (mag(signedDistance[cellLabel]) + mag(signedDistance[interpolationCellLabel]));
        }
    }
}

// TODO: this probably needs some serious optimization (TT)
std::vector<label> frontCompactDivGradCurvatureModel::findNeighbourCells(
    const fvMesh& mesh,
    const label& cellLabel
) const
{
    std::vector<label> neighbourCells{};

    const auto& vertexLabels = mesh.cellPoints()[cellLabel];
    const auto& vertexToCell = mesh.pointCells();

    forAll(vertexLabels, I)
    {
        const auto& connectedCells = vertexToCell[vertexLabels[I]];

        forAll(connectedCells, K)
        {
            if (connectedCells[K] != cellLabel)
            {
                neighbourCells.push_back(connectedCells[K]);
            }
        }
    }

    // Remove duplicate entries
    std::sort(neighbourCells.begin(), neighbourCells.end());
    
    std::vector<label>::iterator newEnd = 
                std::unique(neighbourCells.begin(), neighbourCells.end());

    neighbourCells.resize(std::distance(neighbourCells.begin(), newEnd));

    return neighbourCells;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontCompactDivGradCurvatureModel::frontCompactDivGradCurvatureModel(const dictionary& configDict)
    :
        frontCurvatureModel(configDict),
        distanceCorrection_{configDict.lookup("distanceCorrection")}
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

tmp<volScalarField> frontCompactDivGradCurvatureModel::cellCurvature(
    const fvMesh& mesh, 
    const triSurfaceFront& frontMesh
) const
{
    const Time& runTime = mesh.time();  

    triSurfaceScalarField curvatureBuffer
    {
        IOobject(
            "curvatureBuffer",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        frontMesh,
        dimensionedScalar{
            "zero",
            dimless/dimLength,
            0
        }
    };

    const auto& communication = mesh.lookupObject<lentCommunication>(
                                    lentCommunication::registeredName(frontMesh, mesh)
                                ); 
    const auto& triangleToCell = communication.triangleToCell();
    auto cellCurvatureFieldTmp = frontCurvatureModel::cellCurvature(mesh, frontMesh);
    auto& cellCurvatureField = cellCurvatureFieldTmp.ref();


    const volScalarField& curvatureInputField = 
        mesh.lookupObject<volScalarField>(curvatureInputFieldName()); 

    // Apply corrections
    if (distanceCorrection_ == "sphere")
    {
        applySphereCorrection(mesh, frontMesh, curvatureInputField, cellCurvatureField);
    }
    else if (distanceCorrection_ == "interpolation")
    {
        applyInterpolationCorrection(mesh, frontMesh, curvatureInputField, cellCurvatureField);
    }
    else if (distanceCorrection_ == "off")
    {
    }
    else
    {
        FatalErrorIn (
            "frontCompactDivGradCurvatureModel::cellCurvature(front, mesh)"
        )   << "Unknown distance correction "
            << distanceCorrection_
            << " for the keyword 'distanceCorrection'.\n"
            << "Valid options are 'sphere', 'interpolation' or 'off'."
            << exit(FatalError);
    }

    // Write to temporary front field to communicate curvature back to
    // Eulerian mesh
    forAll(triangleToCell, index)
    {
        curvatureBuffer[index] = cellCurvatureField[triangleToCell[index]];
    }

    // Propagate compact curvature to narrow band cells
    cellCurvatureField *= 0.0;

    const auto& cellToTriangle = communication.cellsTriangleNearest();

    forAll(cellToTriangle, I)
    {
        const auto& hitObject = cellToTriangle[I];
        
        if (hitObject.hit())
        {
            cellCurvatureField[I] = curvatureBuffer[hitObject.index()];
        }
    }

    return cellCurvatureFieldTmp;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
