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
    Foam::lentInterpolationTemplates

SourceFiles
    lentInterpolationTemplates.H

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
        Front / Mesh communication maps.  

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationCellPoint.H"

namespace Foam {
namespace FrontTracking {

template<typename Type>  
void lentInterpolation::interpolate(
    const DimensionedField<Type, triSurfaceFrontPointGeoMesh>& vertexField,
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
    tmp<DimensionedField<Type, triSurfaceFrontPointGeoMesh>> vertexFieldTmp,
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

    interpolationCellPoint<Type> barycentric(cellField);
    const auto& vertexToCell = communication.vertexToCell(); 
    const pointField& vertices = front.points();

    forAll(vertexToCell, vertexI)
    {
        vertexField[vertexI] = barycentric.interpolate(
            vertices[vertexI],
            vertexToCell[vertexI] 
        );
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
