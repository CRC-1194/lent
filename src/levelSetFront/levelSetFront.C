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
    tomislav.maric@gmx.com

\*---------------------------------------------------------------------------*/

#include "levelSetFront.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * //

namespace Foam
{
    namespace frontTracking
    {
        defineTypeNameAndDebug(levelSetFront, 0);
    }
}

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * //

namespace Foam {
    namespace frontTracking {

//void levelSetFront::computeIsoSurface(
    //const volScalarField& cellsToElementsDist, 
    //const scalarField& pointsToElementsDist, 
    //const bool regularise, 
    //const scalar mergeTol
//)
//{
    //*this = isoSurface 
    //(
        //cellsToElementsDist, 
        //pointsToElementsDist, 
        //0, 
        //regularise,
        //mergeTol
    //);
//}

//void levelSetFront::write(label index)
//{
    //// Separate the file name and the extension.
    //fileName baseName = file_.name(true);

    //string indexString = Foam::name(index);

    //// Pad the base name with zeros
    //std::string paddedZeros = std::string (
        //prependZeros_ - indexString.size(), 
        //'0'
    //);

    //// Append the index string to the padded name.
    //paddedZeros.append(indexString);

    //// Write the front in the instance directory under the new name.
    
    //// TODO: generalize the IO for file formats 
    //fileName finalName = directory_ + "/" + 
        ////baseName + "-" + paddedZeros + "." + name_.ext(); 
        //baseName + "-" + paddedZeros + ".vtk"; 

    //write(finalName);
//}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//levelSetFront::levelSetFront(
    //const IOobject& io, 
    //label prependZeros
//)
//:
    //moving_(false), 
    //changing_(false), 
    //name_(io.name()), 
    //instance_(io.instance()), 
    //prependZeros_(prependZeros), 
    //frontPtr_()
//{
//}

levelSetFront::levelSetFront(
    const IOobject& io, 
    label prependZeros
)
:
    regIOobject(io), 
    triSurface(io.filePath())
{
}


//levelSetFront::levelSetFront(const volScalarField& psi,
                                                  //const scalarField& psiPoint)
    //: 
        //triSurfaceMesh(), 
        //moving_(false), 
        //changing_(false)
//{
    //// Reconstruct the interface as a level 0 of the distance field.
    //computeIsoSurface(psi, psiPoint); 
//}


//levelSetFront::levelSetFront(const levelSetFront& rhs)
//:
    //triSurfaceMesh(rhs),
    //moving_(rhs.moving_),
    //changing_(rhs.changing_)
//{

//}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//levelSetFront::~levelSetFront()
//{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//bool levelSetFront::isMoving() const
//{
    //return moving_;
//}

//void levelSetFront::setMoving(bool b)
//{
    //moving_ = b;
//}

//bool levelSetFront::isChanging() const
//{
    //return changing_;
//}
//void levelSetFront::setChanging(bool b)
//{
    //changing_ = b;
//}

//void levelSetFront::reconstruct(
    //const volScalarField& cellsToElementsDist, 
    //const scalarField& pointsToElementsDist, 
    //const bool regularise, 
    //const scalar mergeTol
//)
//{
    //computeIsoSurface (
        //cellsToElementsDist, 
        //pointsToElementsDist, 
        //regularise, mergeTol
    //);
//}

//void levelSetFront::reconstruct(
    //const volScalarField& cellsToElementsDist, 
    //const bool regularise, 
    //const scalar mergeTol
//)
//{
    //const fvMesh& mesh = cellsToElementsDist.mesh(); 
    //volPointInterpolation pInter (mesh);

    //tmp<pointScalarField> pointsToElementsDistTmp = pInter.interpolate (
        //cellsToElementsDist
    //);

    //const scalarField& pointsToElementsDist = pointsToElementsDistTmp();

    //computeIsoSurface (
        //cellsToElementsDist, 
        //pointsToElementsDist, 
        //regularise, mergeTol
    //);
//}

//void levelSetFront::move(vector deltaV)
//{
    //executeMovePoints(deltaV);
//}

//void levelSetFront::move(const vectorField& deltaV)
//{
    //executeMovePoints(deltaV);
//}

//void levelSetFront::write(const Time& runTime)
//{
    //// If the time is the output time.
    //if (runTime.outputTime())
    //{
        //write(runTime.timeIndex()); 
    //}
//}

//void levelSetFront::writeNow(const Time& runTime)
//{
    //// If the time is the output time.
    //write(runTime.timeIndex()); 
//}

bool levelSetFront::write() const
{
    Info << "write() called" << endl;

    return true;
}

bool levelSetFront::writeData(Foam::Ostream& os) const
{
    Info << "writeData() called" << endl;
    return true; 
}



// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //

//void levelSetFront::operator=(const levelSetFront& rhs)
//{
    //if (&rhs != this)
    //{
        //// TODO: add code
    //}
//}

//void levelSetFront::operator=(const isoSurface& rhs)
//{
    //*this = static_cast<const triSurface&> (rhs); 
//}

// ************************************************************************* //

} // End namespace frontTracking

} // End namespace Foam
