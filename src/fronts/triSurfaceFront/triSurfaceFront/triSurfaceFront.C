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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceFront::triSurfaceFront(
    const IOobject& io,
    word writeFormat,
    label prependZeros
)
:
    regIOobject(io),
    triSurface(),
    writeFormat_(writeFormat),
    prependZeros_(prependZeros)
{
    fileName file = existingFrontFileName(io);

    static_cast<triSurface&>(*this) = triSurface(file);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

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
