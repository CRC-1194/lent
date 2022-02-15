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
    Foam::triSurfaceFront

SourceFiles
    triSurfaceFront.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    A registered triSurface that can be used with tmp<>.

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


#include "triSurfaceFront.H"
#include "volPointInterpolation.H"
#include "fvcGrad.H"

#include <algorithm>

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

void triSurfaceFront::computeTwoRingNeighbours()
{
    twoRingNeighbours_.resize(this->localPoints().size());

    const auto& pointToEdges = this->pointEdges();

    forAll(pointToEdges, I)
    {
        std::vector<label> neighbours = neighbourPoints(I);

        auto oneRingNeighbours = neighbours;

        for (const auto& pointID : oneRingNeighbours)
        {
            // Exclude the current point I since it is not a neighbour
            // of itself
            auto additionalNeighbours = neighbourPoints(pointID, I);

            for (auto neighbourID : additionalNeighbours)
            {
                neighbours.push_back(neighbourID);
            }
        }

        std::sort(neighbours.begin(), neighbours.end());

        // Remove duplicates
        auto endIt = std::unique(neighbours.begin(), neighbours.end());
        neighbours.resize(std::distance(neighbours.begin(), endIt));

        twoRingNeighbours_[I] = neighbours;
    }
}

std::vector<label> triSurfaceFront::neighbourPoints(const label& pointLabel, const label& exclude) const
{
    std::vector<label> neighbourLabels{};

    const auto& connectedEdges = this->pointEdges()[pointLabel];
    const auto& edges = this->edges();

    forAll(connectedEdges, I)
    {
        const auto& anEdge = edges[connectedEdges[I]];

        forAll(anEdge, K)
        {
            if (anEdge[K] != pointLabel && anEdge[K] != exclude)
            {
                neighbourLabels.push_back(anEdge[K]);
            }
        }
    }

    return neighbourLabels;
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
    prependZeros_(prependZeros),
    twoRingNeighbours_{}
{
    // FIXME: Output requires special handling for a parallel run. TM.
    // Create the "front" directory in the case sub-directory. 
    mkDir
    (
        IOobject::rootPath() + "/" + 
        IOobject::caseName() + "/" 
        + instance()
    ); 

    // FIXME: This is ugly, clean this up. TM.
    fileName initialFileName = IOobject::rootPath() + "/" + 
        IOobject::caseName() + "/" +
        this->time().constant() + "/" + 
        IOobject::name() + "." + readFormat_;  

    // Construct the triSurface from the current file. 
    static_cast<triSurface&>(*this) = triSurface(initialFileName);

    computeTwoRingNeighbours();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool triSurfaceFront::write(const bool) const
{
    if (writeFormat_ == "vtk")
    {
        OFstream output(actualFileName());
        writeVTKWithFields(output);
    }
    else
    {
        triSurface::write(actualFileName());
    }

    return true;
}

bool triSurfaceFront::writeData(Foam::Ostream& os) const
{
    if (writeFormat_ == "vtk")
    {
        writeVTKWithFields(os);
    }
    else
    {
        triSurface::write(os);
    }

    return true;
}

bool triSurfaceFront::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType,
    const bool
) const
{
    if (writeFormat_ == "vtk")
    {
        OFstream output(actualFileName());
        writeVTKWithFields(output);
    }
    else
    {
        triSurface::write(actualFileName());
    }

    return true;
}

void triSurfaceFront::displace(const Field<Vector<double> >& displacements)
{
    auto& frontVertices = this->storedPoints();
    frontVertices += displacements;

    // This deletes all demand driven geometrical data (e.g. face normals).
    // Thus, it is recomputed once it is called for. However, the memory
    // is also deallocated.
    // TODO: When profiling Lent, take a look at the performance impact of
    // the continued allocation / deallocation
    clearGeom();
}

point triSurfaceFront::geometricCentre() const
{
    point geoCentre{0,0,0};

    for (const auto& vertex : this->points())
    {
        geoCentre += vertex;
    }

    geoCentre /= this->points().size();

    return geoCentre;
}

scalar triSurfaceFront::convexFrontVolume() const
{
    // NOTE: as the name suggests, this simple approach only works
    // for convex interfaces (TT)
    // Just uses a tetrahedral decompostion of the enclosed volume and
    // the triple product
    scalar volume{0.0};

    const auto& V = this->points();
    const List<labelledTri>& triangles = static_cast<const List<labelledTri>& > (*this);
    const auto gc = geometricCentre();

    for (const auto& T : triangles)
    {
        volume += mag(((V[T[1]] - V[T[0]]) ^ (V[T[2]] - V[T[0]])) & (gc - V[T[0]]));
    }

    volume /= 6.0;

    return volume;
}

// * * * * * * * * * * * * * * Member Operators * * * * * * * * * * * * * * //
void triSurfaceFront::operator=(const triSurface& rhs)
{
    triSurface::operator=(rhs);  

    computeTwoRingNeighbours();
}

void triSurfaceFront::operator=(const isoSurfacePoint& rhs)
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

    computeTwoRingNeighbours();
}

// ************************************************************************* //

} // End namespace FrontTracking

} // End namespace Foam
