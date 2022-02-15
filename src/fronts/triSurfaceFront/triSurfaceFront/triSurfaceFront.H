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


#ifndef triSurfaceFront_H
#define triSurfaceFront_H

#include "isoSurfacePoint.H"
#include "triSurface.H"
#include "error.H"

#include <vector>

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

    // Topological information required for parabolda fitting
    std::vector<std::vector<label>> twoRingNeighbours_;

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

    // Write fields registered to the front to the vtk file to enable visualization
    // in paraview
    void writeVTKWithFields(Ostream& os) const;

    void computeTwoRingNeighbours();
    std::vector<label> neighbourPoints(const label& pointLabel, const label& exclude=-1) const;

public:

    TypeName ("triSurfaceFront");

    // Constructors
    triSurfaceFront(const IOobject& io, word readFormat="stl", word writeFormat="vtk", label prependZeros=8);

    //- Destructor
    virtual ~triSurfaceFront() {};

    // Member Functions
    bool write(const bool valid = true) const override;

    label nFaces() const 
    {
        return localFaces().size(); 
    }

    void displace(const Field<Vector<double>>& displacements);

    bool writeObject
    (
        IOstream::streamFormat fmt,
        IOstream::versionNumber ver,
        IOstream::compressionType cmp,
        const bool valid = true
    ) const override;

    bool writeData(Foam::Ostream& os) const override;

    const std::vector<std::vector<label>>& twoRingNeighbours() const
    {
        return twoRingNeighbours_;
    }

    point geometricCentre() const;
    scalar convexFrontVolume() const;

    // Member operators
    void operator=(const triSurfaceFront& rhs);
    void operator=(const triSurface& rhs);

    void operator=(const isoSurfacePoint& rhs);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
