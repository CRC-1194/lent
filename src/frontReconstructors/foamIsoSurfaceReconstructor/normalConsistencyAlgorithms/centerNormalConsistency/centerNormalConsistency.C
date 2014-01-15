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

Authors
    Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    Mathematical Modeling and Analysis
    Center of Smart Interfaces, TU Darmstadt

\*---------------------------------------------------------------------------*/

#include "centerNormalConsistency.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H" 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking { 

    defineTypeNameAndDebug(centerNormalConsistency, 0); 
    addToRunTimeSelectionTable(normalConsistency, centerNormalConsistency, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

centerNormalConsistency::centerNormalConsistency(const dictionary& configDict)
    :
        normalConsistency(configDict)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void centerNormalConsistency::makeFrontNormalsConsistent(
    triSurface& front, 
    const labelList& elementCells, 
    const volScalarField& signedDistance
) const
{
    // Gradient based centerNormal consistency algorithm.
    volVectorField distGrad = fvc::grad(signedDistance); 

    List<labelledTri>& elements = static_cast<List<labelledTri>& > (front);  

    const vectorField& elementNormals = front.faceNormals();  

    // For all faces 
    forAll (elements, E)
    {
        scalar gradMag = mag(distGrad[elementCells[E]]); 

        if (gradMag >= SMALL)
        {
            distGrad[elementCells[E]] /= gradMag; 
            
            scalar centerNormalMag = mag(elementNormals[E]); 

            if (centerNormalMag > SMALL)
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
