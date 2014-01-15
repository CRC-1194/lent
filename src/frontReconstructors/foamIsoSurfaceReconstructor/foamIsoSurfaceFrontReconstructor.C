/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "foamIsoSurfaceFrontReconstructor.H"
#include "addToRunTimeSelectionTable.H"
#include "isoSurface.H"
#include "fvcGrad.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(foamIsoSurfaceFrontReconstructor, 0); 
    addToRunTimeSelectionTable(foamIsoSurfaceFrontReconstructor, foamIsoSurfaceFrontReconstructor, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

foamIsoSurfaceFrontReconstructor::foamIsoSurfaceFrontReconstructor(
   const dictionary& configDict
)
:
    frontReconstructor(configDict), 
    mergeTolerance_(readScalar(configDict.lookup("mergeTolerance"))), 
    regularize_(configDict.lookup("regularization")),
    consistencyAlgTmp_(normalConsistency::New(configDict.subDict("normalConsistency"))) 
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

foamIsoSurfaceFrontReconstructor::~foamIsoSurfaceFrontReconstructor()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

labelList foamIsoSurfaceFrontReconstructor::reconstructFront(
    triSurfaceFront& front, 
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance
) const
{
    isoSurface iso (
        signedDistance, 
        pointSignedDistance, 
        0, 
        regularize_,
        mergeTolerance_
    );

    consistencyAlgTmp_->makeFrontNormalsConsistent(
        iso, 
        iso.meshCells(), 
        signedDistance
    );

    front = iso; 

    return iso.meshCells(); 
} 
void foamIsoSurfaceFrontReconstructor::forceConsistentNormalOrientation(
    isoSurface& iso, 
    const volScalarField& signedDistance
) const
{
    volVectorField distGrad = fvc::grad(signedDistance); 

    // Get non-const access to elements.
    List<labelledTri>& elements = static_cast<List<labelledTri>& > (iso);  

    // Get the cells. 
    const labelList& elementCells = iso.meshCells(); 

    // Get the element normals. 
    const vectorField& elementNormals = iso.faceNormals();  

    // For all faces 
    forAll (elements, E)
    {
        // Normalize the distance gradient to get only the direction. 
        scalar gradMag = mag(distGrad[elementCells[E]]); 

        if (gradMag >= SMALL)
        {
            distGrad[elementCells[E]] /= gradMag; 
            
            scalar normalMag = mag(elementNormals[E]); 

            if (normalMag > SMALL)
            {
                vector elementNormal = elementNormals[E] / mag(elementNormals[E]); 

                if ((elementNormal & distGrad[elementCells[E]]) < 0)
                {
                    elements[E].flip(); 
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
