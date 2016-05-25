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
    Foam::tetFillingLevelAlternativeMarkerFieldModel

SourceFiles
    tetFillingLevelAlternativeMarkerFieldModel.C

Author
    Tobias Tolle    tolle@csi.tu-darmstadt.de

Description
    Computes the markerfield from the vertex based and cell centre based 
    signed distances by decomposition of a cell into tetrahedrons and calculating
    the filling level for each tetrahedron as described in

        "From level set to volume of fluid and back again at second-order 
         accuracy" Detrixhe and Aslam

    This class uses an alternative approach for the triangulation of the faces
    which does not introduce additional points and thus avoids interpolation
    errors regardign the signed distance.

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

#ifndef tetFillingLevelAlternativeMarkerFieldModel_H
#define tetFillingLevelAlternativeMarkerFieldModel_H

#include "typeInfo.H"

#include "tetFillingLevelMarkerFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
         Class tetFillingLevelAlternativeMarkerFieldModel Declaration
\*---------------------------------------------------------------------------*/

class tetFillingLevelAlternativeMarkerFieldModel
:
    public tetFillingLevelMarkerFieldModel
{

public:

    TypeName ("tetFillingLevelAlternative");

    // Constructors
    using tetFillingLevelMarkerFieldModel::tetFillingLevelMarkerFieldModel;

    // Member Functions
    virtual List<tetrahedron> tetDecomposition(const point& cellCentre,
                                               const face& cellFace,
                                               const pointField& points,
                                               const pointScalarField& pointDistance,
                                               scalar centreDistance) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
