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

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

#include "levelSetFront.H"
#include "volPointInterpolation.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * //

namespace Foam
{
    namespace frontTracking
    {
        defineTypeNameAndDebug(levelSetFront, 0);
    }
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

void Foam::frontTracking::levelSetFront::computeIsoSurface (
    const volScalarField& cellsToElementsDist, 
    const scalarField& pointsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    *this = isoSurface (
        cellsToElementsDist, 
        pointsToElementsDist, 
        0, 
        regularise,
        mergeTol
    );
}

void Foam::frontTracking::levelSetFront::write(label index)
{
    // Separate the file name and the extension.
    fileName baseName = name_.name(true);

    string indexString = Foam::name(index);

    // Pad the base name with zeros
    std::string paddedZeros = std::string (
        prependZeros_ - indexString.size(), 
        '0'
    );
        //std::string dest = std::string( 10 - original_string.size() , '0').append( original_string);

    // Append the index string to the padded name.
    paddedZeros.append(indexString);

    // Write the front in the instance directory under the new name.
    
    // FIXME: this will run on POSIX only. Implement system separators.
    fileName finalName = instance_ + "/" + 
        baseName + "-" + paddedZeros + "." + name_.ext(); 

    triSurface::write(finalName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::frontTracking::levelSetFront::levelSetFront (
    const IOobject& io, 
    label prependZeros
)
:
    triSurfaceMesh(io), 
    moving_(false), 
    changing_(false), 
    name_(io.name()), 
    instance_(io.instance()), 
    prependZeros_(prependZeros)
{
    Pout << "levelSetFront::levelSetFront(const IOobject& io)" << endl;
}

//Foam::frontTracking::levelSetFront::levelSetFront(const volScalarField& psi,
                                                  //const scalarField& psiPoint)
    //: 
        //triSurfaceMesh(), 
        //moving_(false), 
        //changing_(false)
//{
    //// Reconstruct the interface as a level 0 of the distance field.
    //computeIsoSurface(psi, psiPoint); 
//}


//Foam::frontTracking::levelSetFront::levelSetFront(const levelSetFront& rhs)
//:
    //triSurfaceMesh(rhs),
    //moving_(rhs.moving_),
    //changing_(rhs.changing_)
//{

//}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//Foam::frontTracking::levelSetFront::~levelSetFront()
//{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::frontTracking::levelSetFront::isMoving() const
{
    return moving_;
}

void Foam::frontTracking::levelSetFront::setMoving(bool b)
{
    moving_ = b;
}

bool Foam::frontTracking::levelSetFront::isChanging() const
{
    return changing_;
}
void Foam::frontTracking::levelSetFront::setChanging(bool b)
{
    changing_ = b;
}

void Foam::frontTracking::levelSetFront::reconstruct (
    const volScalarField& cellsToElementsDist, 
    const scalarField& pointsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    Pout << "Foam::frontTracking::levelSetFront::reconstruct (\n"
         << "    const volScalarField& cellsToElementsDist, \n"
         << "    const scalarField& pointsToElementsDist, \n"
         << "    const bool regularise, \n"
         << "    const scalar mergeTol\n)" << endl;

    //computeIsoSurface (
        //cellsToElementsDist, 
        //pointsToElementsDist, 
        //regularise, mergeTol
    //);

    //setChanging(true);
}

void Foam::frontTracking::levelSetFront::reconstruct (
    const volScalarField& cellsToElementsDist, 
    const bool regularise, 
    const scalar mergeTol
)
{
    Pout << "Foam::frontTracking::levelSetFront::reconstruct (\n"
         << "    const volScalarField& cellsToElementsDist, \n"
         << "    const bool regularise, \n"
         << "    const scalar mergeTol\n)" << endl;

    //const fvMesh& mesh = cellsToElementsDist.mesh(); 
    //volPointInterpolation pInter (mesh);

    //tmp<pointScalarField> pointsToElementsDistTmp = pInter.interpolate (
        //cellsToElementsDist
    //);

    //const pointScalarField& pointsToElementsDist = pointsToElementsDistTmp();

    //computeIsoSurface (
        //cellsToElementsDist, 
        //pointsToElementsDist, 
        //regularise, mergeTol
    //);

    //setChanging(true);
}

void Foam::frontTracking::levelSetFront::move(const vectorField& Dv)
{
    // For all front points
        // Add displacement to the point vector

    setMoving(true);
}

void Foam::frontTracking::levelSetFront::write(const Time& runTime)
{

    // Get the control dictionary.
    // Get the write option from the control dictionary.
    
    // If write option is set to time step 
        // If the time index is divisible by the write interval 
        
            // Write the front in the instance directory with the index 
            // coming from the time-step index.
            write(runTime.timeIndex()); 

    // Else If the write option is set to time value 
        // If the time value is right (check how Time does this)
        
            // Write the fron int he instance directory with the index
            // coming from the time step index.

}



// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

void Foam::frontTracking::levelSetFront::operator=(const isoSurface& rhs)
{
    triSurface::operator=(rhs);
}

// ************************************************************************* //
