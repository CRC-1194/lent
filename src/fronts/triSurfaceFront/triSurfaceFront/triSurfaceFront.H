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
    triSurfaceFront.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    A registered triSurface that can be used with tmp<>.

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


#ifndef triSurfaceFront_H
#define triSurfaceFront_H

#include "isoSurface.H"
#include "error.H"

#include "frontConstructor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

/*---------------------------------------------------------------------------*\
                         Class triSurfaceFront Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceFront
    :
        public refCount,
        public objectRegistry, 
        //public regIOobject, 
        public triSurface
{
    // Read format.
    word readFormat_;

    // Write format.
    word writeFormat_;

    // Prepend zeros.
    label prependZeros_;

    // Compute the iso-surface.
    void computeIsoSurface(
        const volScalarField& cellsToElementsDist,
        const scalarField& pointsToElementsDist,
        const bool regularise=false,
        const scalar mergeTol=1e-06
    );

    // Create the file name with IOobject::name, extension and timeIndex
    // paddeed with zeros. 
    fileName zeroPaddedFileName(word extension) const;

    template<class Displacement>
    void executeMovePoints(const Displacement& d);

    // Get the actual file name of the stored front. IOobject::name + time
    // index prepended with zeros + extension. 
    fileName actualFileName() const;

public:

    TypeName ("triSurfaceFront");

    // Constructors
    triSurfaceFront(const IOobject& io, word readFormat="stl", word writeFormat="vtk", label prependZeros=8);

    //- Destructor
    virtual ~triSurfaceFront() {};

    // Member Functions
    bool write() const;

    label nFaces() const 
    {
        return localFaces().size(); 
    }

    // Return non-const access to face normals.
    pointField& storedFaceNormals()
    {
        return const_cast<pointField&>(faceNormals());
    } 

    bool writeObject
    (
        IOstream::streamFormat fmt,
        IOstream::versionNumber ver,
        IOstream::compressionType cmp
    ) const;

    bool writeData(Foam::Ostream& os) const;

    // Member operators
    void operator=(const triSurfaceFront& rhs);
    void operator=(const isoSurface& rhs);
    void operator=(const triSurface& rhs);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
