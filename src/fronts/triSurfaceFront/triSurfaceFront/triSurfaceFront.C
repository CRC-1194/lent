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
    Foam::triSurfaceFront

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

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

fileName triSurfaceFront::zeroPaddedFileName(word writeFormat) const
{
    // Get the time index from Time.
    string index = Foam::name(IOobject::time().timeIndex());

    // Pad the index string with zeros.
    std::string paddedIndex = std::string(prependZeros_ - index.size(), '0');

    // Append the index string to the padded name.
    paddedIndex.append(index);

    // Create the final name from the path, the IOobject name and the
    // extension.
    return (
        IOobject::path() + "/" + IOobject::name() +  "-" 
        + paddedIndex + "." + writeFormat
    );
}

fileName triSurfaceFront::actualFileName() const
{
    // Get the initial full path name from the IOobject. 
    fileName actualFileName = IOobject::path() + 
        "/" + IOobject::name() + "." + writeFormat_; 

    // Padd the name with zeros and add extension to the IOobject file name
    // so that ParaView can open the temporal file sequence.
    actualFileName = zeroPaddedFileName(writeFormat_);

    return actualFileName;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceFront::triSurfaceFront(
    const IOobject& io,
    word readFormat,
    word writeFormat,
    label prependZeros
)
:
    objectRegistry(io),
    //regIOobject(io), 
    triSurface(),
    readFormat_(readFormat),
    writeFormat_(writeFormat),
    prependZeros_(prependZeros)
{

    // FIXME: Work here to re-start the computation from latestTime.  Get the
    // current file name of the front from the IOobject using runTime and
    // readFormat. . TM.  
    fileName initialFileName = IOobject::path() + 
        "/" + IOobject::name() + "." + readFormat_; 

    // Construct the triSurface from the current file. 
    static_cast<triSurface&>(*this) = triSurface(initialFileName);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool triSurfaceFront::write() const
{
    triSurface::write(actualFileName());

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
    triSurface::write(actualFileName());

    return true;
}

// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //
void triSurfaceFront::operator=(const triSurface& rhs)
{
    triSurface::operator=(rhs);  
}

void triSurfaceFront::operator=(const isoSurface& rhs)
{
    this->clearOut(); 

    auto& thisPoints = this->storedPoints(); 
    auto& thisFaces = this->storedFaces(); 

    const auto& rhsPoints = rhs.localPoints();  
    const auto& rhsFaces = rhs.localFaces(); 

    thisPoints = rhsPoints; 
    thisFaces.resize(rhsFaces.size());
    forAll(rhsFaces, faceI)
    {
        auto& thisFace = thisFaces[faceI]; 
        const auto& rhsFace = rhsFaces[faceI]; 
        thisFace.resize(rhsFace.size()); 
        forAll(rhsFace, pointI)
        {
            thisFace[pointI] = rhsFace[pointI]; 
            thisFace.region() = 0;  
        }
    }
}

// ************************************************************************* //

} // End namespace FrontTracking

} // End namespace Foam
