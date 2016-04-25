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
    Foam::tetFillingLevelMarkerFieldModel

SourceFiles
    tetFillingLevelMarkerFieldModel.C

Author
    Tobias Tolle    tolle@csi.tu-darmstadt.de

Description
    Computes the markerfield from the vertex based and cell centre based 
    signed distances by decomposition of a cell into tetrahedrons and calculating
    the filling level for each tetrahedron as described in

        "From level set to volume of fluid and back again at second-order 
         accuracy" Detrixhe and Aslam

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

#ifndef tetFillingLevelMarkerFieldModel_H
#define tetFillingLevelMarkerFieldModel_H

#include "typeInfo.H"
#include "pointFields.H"
#include "volFieldsFwd.H"

#include "markerFieldModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class tetFillingLevelMarkerFieldModel Declaration
\*---------------------------------------------------------------------------*/

class tetFillingLevelMarkerFieldModel
:
    public markerFieldModel
{
    // Private data
    const word pointDistFieldName_;


    // Private Member Functions
    void tagInterfaceCells(volScalarField& markerField) const;
    void setBulkMarkerField(volScalarField& markerField) const;
    scalar fillingLevel(const label& cellID, const face& cellFace,
                        const fvMesh& mesh) const;
    scalar tetrahedralVolume(const point& pointA, const point& pointB,
                             const point& pointC, const point& pointD) const;
    scalar volumeFraction(List<scalar>& d) const;
    label sortDistances(List<scalar>& distances) const;


public:

    TypeName ("tetFillingLevel");

    // Constructors
    tetFillingLevelMarkerFieldModel(const dictionary& configDict);


    // Member Functions
    virtual void calcMarkerField(volScalarField& markerField) const;

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
