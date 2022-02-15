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
    Foam::tetFillingLevelMarkerFieldModel

SourceFiles
    tetFillingLevelMarkerFieldModel.C

Author
    Tobias Tolle (tolle@mma.tu-darmstadt.de)

Description
    Computes the markerfield from the vertex based and cell centre based 
    signed distances by decomposition of a cell into tetrahedrons and calculating
    the filling level for each tetrahedron as described in

        "From level set to volume of fluid and back again at second-order 
         accuracy" Detrixhe and Aslam

    This class uses a barycentric decomposition of the cell faces into triangles.

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

#ifndef tetFillingLevelMarkerFieldModel_H
#define tetFillingLevelMarkerFieldModel_H

#include "typeInfo.H"
#include "pointFields.H"
#include "SortableList.H"
#include "volFieldsFwd.H"

#include "markerFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                Class tetFillingLevelMarkerFieldModel Declaration
\*---------------------------------------------------------------------------*/

struct tetrahedron
{
    FixedList<point, 4> vertex;
    SortableList<scalar> distance;
};

class tetFillingLevelMarkerFieldModel
:
    public markerFieldModel
{
    // Private data
    const word pointDistFieldName_;
    const label nSmoothingSteps_;

public:

    TypeName ("tetFillingLevel");

    // Constructors
    tetFillingLevelMarkerFieldModel(const dictionary& configDict);

    // Member Functions
    virtual void calcMarkerField(volScalarField& markerField) const;

    void tagInterfaceCells(volScalarField& markerField) const;
    void setBulkMarkerField(volScalarField& markerField) const;
    virtual List<tetrahedron> tetDecomposition(const point& cellCentre,
                                               const face& cellFace,
                                               const pointField& points,
                                               const pointScalarField& pointDistance,
                                               scalar centreDistance) const;
    scalar fillingLevel(const label& cellID, const face& cellFace,
                         const fvMesh& mesh) const;
    scalar tetrahedralVolume(const tetrahedron& tet) const;
    scalar volumeFraction(SortableList<scalar>& d) const;
    label countNegativeEntries(List<scalar>& distance) const;
    void smoothMarkerField(volScalarField& markerField) const;
    void cutOverUnderShoots(volScalarField& markerField) const;

    word pointDistFieldName() const
    {
        return pointDistFieldName_; 
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
