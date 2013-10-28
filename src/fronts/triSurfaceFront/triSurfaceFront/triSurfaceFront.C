/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Tomislav Maric
    maric<<at>>csi<<dot>>tu<<minus>>darmstadt<<dot>>de
    tomislav<<dot>>maric<<at>>gmx<<dot>>com

\*---------------------------------------------------------------------------*/

#include "triSurfaceFront.H"
#include "volPointInterpolation.H"
#include "fvcGrad.H" 

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * //

namespace Foam
{
    namespace FrontTracking
    {
        defineTypeNameAndDebug(triSurfaceFront, 0);
    }
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

namespace Foam {
    namespace FrontTracking {

void triSurfaceFront::computeIsoSurface(
    const volScalarField& cellsToElementsDist, 
    const scalarField& pointsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    isoSurface reconstructedFront = isoSurface 
    (
        cellsToElementsDist, 
        pointsToElementsDist, 
        0, 
        regularise,
        mergeTol
    );

    meshCells_ = reconstructedFront.meshCells(); 
    
    // Copy the data to the front. 
    *this = reconstructedFront; 

    // FIXME: careful, BUG inside: normals are not consistent 
    forceConsistentNormalOrientation(cellsToElementsDist);  
}

fileName triSurfaceFront::zeroPaddedFileName(word extension) const
{
    // Separate the file name and the extension.
    fileName file = IOobject::name();
    fileName baseName = file.name(true); 

    string indexString = Foam::name(IOobject::time().timeIndex());
    // Pad the base name with zeros
    std::string paddedZeros = std::string (
        prependZeros_ - indexString.size(), 
        '0'
    );

    // Append the index string to the padded name.
    paddedZeros.append(indexString);

    fileName finalName = path() + "/" + baseName + "-" +
        paddedZeros + "." + extension; 

    return finalName; 
}

fileName triSurfaceFront::existingFrontFileName(const IOobject& io)
{
    fileName modifiedFileName = io.filePath();
    
    const Time& runTime = io.time(); 

    if (runTime.timeIndex() > 0)
    {
        modifiedFileName = zeroPaddedFileName(writeFormat_);  
    }

    return modifiedFileName;
} 

void triSurfaceFront::forceConsistentNormalOrientation
(
    //const isoSurface& reconstructedFront, 
    const volScalarField& cellsToElementsDist
)
{
    volVectorField distGrad = fvc::grad(cellsToElementsDist); 

    // Get non-const access to elements.
    List<labelledTri> & elements = storedFaces(); 

    // Get the cells. 
    const labelList& elementCells = meshCells(); 

    // Get the element normals. 
    // FIXME: isoSurface: faceNormals are not updated when an in-place flip of an element.
    const vectorField& elementNormals = faceNormals();  

    // For all faces 
    forAll (elements, E)
    {
        // Normalize the distance gradient to get only the direction. 
        // TODO: remove, debugging only
        scalar gradMag = mag(distGrad[elementCells[E]]); 

        // TODO: remove, debugging only
        if (gradMag >= SMALL)
        {
            distGrad[elementCells[E]] /= gradMag; 
            
            // TODO: remove, debugging only
            scalar normalMag = mag(elementNormals[E]); 

            // TODO: remove, debugging only
            if (normalMag > SMALL)
            {
                vector elementNormal = elementNormals[E] / mag(elementNormals[E]); 

                // If the normal is oriented in the opposite way from the distance gradient. 
                if ((elementNormal & distGrad[elementCells[E]]) < 0)
                {
                    elements[E].flip(); 
                }
            }
        }
    }
};

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceFront::triSurfaceFront(
    const IOobject& io, 
    word writeFormat, 
    label prependZeros
)
:
    regIOobject(io), 
    //triSurface(io.filePath()), 
    triSurface(), 
    meshCells_(), 
    writeFormat_(writeFormat),
    prependZeros_(prependZeros)
{
    fileName file = existingFrontFileName(io);

    static_cast<triSurface&>(*this) = triSurface(file);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void triSurfaceFront::reconstruct(
    const volScalarField& cellsToElementsDist, 
    const scalarField& pointsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    computeIsoSurface (
        cellsToElementsDist, 
        pointsToElementsDist, 
        regularise, mergeTol
    );
}

void triSurfaceFront::reconstruct(
    const volScalarField& cellsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    const fvMesh& mesh = cellsToElementsDist.mesh(); 
    volPointInterpolation pInter (mesh);

    tmp<pointScalarField> pointsToElementsDistTmp = pInter.interpolate (
        cellsToElementsDist
    );

    const scalarField& pointsToElementsDist = pointsToElementsDistTmp();

    computeIsoSurface (
        cellsToElementsDist, 
        pointsToElementsDist, 
        regularise, mergeTol
    );
}

void triSurfaceFront::move(vector deltaV)
{
    executeMovePoints(deltaV);
}

void triSurfaceFront::move(const vectorField& deltaV)
{
    executeMovePoints(deltaV);
}

bool triSurfaceFront::write() const
{
    fileName paddedName = zeroPaddedFileName(writeFormat_);

    triSurface::write(paddedName); 

    return true;
}

bool triSurfaceFront::writeData(Foam::Ostream& os) const
{
    triSurface::write(os); 

    return true; 
}

bool triSurfaceFront::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    triSurface::write(zeroPaddedFileName(writeFormat_)); 

    return true;
}

// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

void triSurfaceFront::operator=(const isoSurface& rhs)
{

    static_cast<triSurface*>(this)->operator=(static_cast<const triSurface&> (rhs)); 
}

// ************************************************************************* //

} // End namespace FrontTracking

} // End namespace Foam
