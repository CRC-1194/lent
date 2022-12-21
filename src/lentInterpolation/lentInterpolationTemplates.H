/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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
    Foam::lentInterpolation

SourceFiles
    lentInterpolation.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    LENT interpolation interface
        
Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationCellPoint.H"
#include "dualKrigingInterpolation.H"
#include "fvcGrad.H"

#include <chrono>

namespace Foam {
namespace FrontTracking {

template<typename Type>  
void lentInterpolation::interpolate(
    const DimensionedField
    <
        Type, 
        triSurfaceFrontPointGeoMesh
    >& vertexField,
    GeometricField<Type, fvPatchField, volMesh>& cellField 
) 
{
    const triSurfaceFront& front = vertexField.mesh(); 
    const fvMesh& mesh = cellField.mesh(); 
    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    // Simply transfer values from the nearest vertex. 
    // FIXME: Implement a barycentrid/volume weighted interpolation here! TM. 
    const auto& vertexToCell = communication.vertexToCell(); 
    forAll(vertexToCell, vertexI)
    {
        cellField[vertexToCell[vertexI]] = vertexField[vertexI]; 
    }
}

template<typename Type>  
void lentInterpolation::interpolate(
    tmp
    <
        DimensionedField<Type, triSurfaceFrontPointGeoMesh>
    > vertexFieldTmp,
    GeometricField<Type, fvPatchField, volMesh>& cellField 
) 
{
    return interpolate(vertexFieldTmp(), cellField); 
}

template<typename Type>  
void lentInterpolation::interpolate(
    const GeometricField<Type, fvPatchField, volMesh>& cellField, 
    DimensionedField<Type, triSurfaceFrontPointGeoMesh>& vertexField
) 
{
    const triSurfaceFront& front = vertexField.mesh(); 
    const fvMesh& mesh = cellField.mesh(); 
    const lentCommunication& communication = 
        mesh.lookupObject<lentCommunication>(
            lentCommunication::registeredName(front,mesh)
    ); 

    if (cellToVertexInterpType_ == "Kriging")
    {
        // Time measurement
        auto start = std::chrono::system_clock::now();

        dualKrigingInterpolation kriging{};

        const auto& vertices = front.points();
        const auto& verticesInCell = communication.interfaceCellToVertices();

        for (const auto& cellToVertices : verticesInCell)
        {
            const auto& cellLabel = cellToVertices.first;
            const auto& vertexLabels = cellToVertices.second;

            kriging.computeDualKrigingParameters(cellLabel, mesh, cellField);

            for (const auto& vLabel : vertexLabels)
            {
                vertexField[vLabel] = kriging.interpolateTo(vertices[vLabel]);
            }

            //kriging.selfTest(cellField, cellLabel);
        }
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsedSeconds = end - start;

        Info << "Computation time Kriging interpolation: " << elapsedSeconds.count() << endl;
    }
    else if (cellToVertexInterpType_ == "extrapolation")
    {
	Info << "Using extrapolation to transfer " << cellField.name() 
		<< " to " << vertexField.name() << endl;
        const auto& vertices = front.points();
        const auto& verticesInCell = communication.interfaceCellToVertices();

	auto cellFieldGradTmp = fvc::grad(cellField, "leastSquares");
	const auto& cellFieldGrad = cellFieldGradTmp();
	const auto& mesh = cellFieldGrad.mesh();
	const auto& cellCenters = mesh.C();

        for (const auto& cellToVertices : verticesInCell)
        {
            const auto& cellLabel = cellToVertices.first;
            const auto& vertexLabels = cellToVertices.second;

            for (const auto& vLabel : vertexLabels)
            {
                vertexField[vLabel] = cellField[cellLabel] + 
			(cellFieldGrad[cellLabel] & (vertices[vLabel] - cellCenters[cellLabel]));
            }
        }
    }
    else
    {
        auto interpolationPtr = 
            interpolation<Type>::New(cellToVertexInterpType_,cellField);  

        const auto& vertexToCell = communication.vertexToCell(); 
        const pointField& vertices = front.points();

        forAll(vertexToCell, vertexI)
        {
            // FIXME: Dynamic dispatch overhead used here.
            // Fix interpolation once we know which type gives best results.
            // TM. 
            vertexField[vertexI] = interpolationPtr->interpolate(
                vertices[vertexI],
                vertexToCell[vertexI] 
            );
        }
    }
}

template<typename Type>  
void lentInterpolation::interpolate(
    tmp<GeometricField<Type, fvPatchField, volMesh>> cellFieldTmp,
    DimensionedField<Type, triSurfaceFrontPointGeoMesh>& vertexField
) 
{
    return interpolate(cellFieldTmp(), vertexField); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
