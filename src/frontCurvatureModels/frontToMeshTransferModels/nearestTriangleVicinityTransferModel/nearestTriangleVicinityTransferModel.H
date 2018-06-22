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
    Foam::nearestTriangleVicinityTransferModel.H

SourceFiles
    nearestTriangleVicinityTransferModel.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description

    Directly transfer the curvature from the front to the cell faces of interface
    cells by finding the closest triangle in a cell neighbourhood.
    The curvature is computed as a weighted average of all triangles contained in
    a ball centred around the barycentre of the closest triangle to the centre of a
    cell face

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

#ifndef nearestTriangleVicinityTransferModel_H
#define nearestTriangleVicinityTransferModel_H

#include "frontToMeshTransferModel.H"

#include "lentCommunication.H"

#include <map>
#include <vector>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class nearestTriangleVicinityTransferModel Declaration
\*---------------------------------------------------------------------------*/

class nearestTriangleVicinityTransferModel
:
    public frontToMeshTransferModel
{
public:

    // Type aliases
    using labelPointPair = std::pair<label, point>;

private:

    // Private data
    const word markerFieldName_;
    scalar searchRadiusCoefficient_;
    mutable scalar searchRadiusSquared_;
    mutable tmp<surfaceScalarField> faceCurvatureTmp_;

    // Buffer for communication data
    mutable std::map<label, std::vector<label>> trianglesInCellNeighbourhood_;

    // Private Member Functions
    void initializeFaceCurvatureField(const fvMesh&) const; 
    void resetCurvature() const;

    void computeTrianglesInCellNeighbourhoodMapping(
                        const fvMesh& mesh,
                        const lentCommunication& communication
                    ) const;

    labelPointPair closestFrontElement(
                        const point& P,
                        const std::vector<label>& triangleIDs,
                        const triSurfaceFront& front
                    ) const;

    label cellContainingClosestElement(
                        const label& triangleID,
                        const lentCommunication& communication
                    ) const;

    std::vector<label> computeTrianglesInVicinity(
                        const point& P,
                        const std::vector<label>& candidates,
                        const triSurfaceFront& front
                    ) const;

    scalar weightedCurvatureAverage(
                        const point& P,
                        const std::vector<label>& triangleIDs,
                        const triSurfaceFrontVectorField& curvatureNormals,
                        const triSurfaceFront& front
                    ) const;

    scalar computeSearchRadiusSquared() const;
    scalar weight(const scalar& distance) const;
    void setSearchRadiusSquared(const fvMesh& mesh) const;
    bool isInterfaceCell(const scalar& volFraction) const;

    // Taken from Dual Kriging interpolation
    void removeDuplicates(std::vector<label>& listOfLabels) const;
    std::vector<label> cellNeighbourhood(const label& cellLabel, const fvMesh& mesh) const;


public:

    // Static data members
    TypeName ("nearestTriangleVicinity");


    // Constructors
    nearestTriangleVicinityTransferModel(const dictionary& configDict);


    //- Destructor
    virtual ~nearestTriangleVicinityTransferModel() = default;


    // Member Functions
    virtual void transferCurvature(
        const triSurfaceFrontVectorField& curvatureNormals,
        const triSurfaceFront& front,
        const fvMesh& mesh
    ) const;

    virtual tmp<volScalarField> cellCurvature() const;  

    virtual tmp<surfaceScalarField> faceCurvature() const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
