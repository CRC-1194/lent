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
    Foam::markerFieldModel

SourceFiles
    markerFieldModel.C

Author
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Abstract base class for the markerField function calculation from a signed
    distance field.
        
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


#ifndef markerFieldModel_H
#define markerFieldModel_H

#include "typeInfo.H"
#include "autoPtr.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class markerFieldModel Declaration
\*---------------------------------------------------------------------------*/

class markerFieldModel
    :
        public refCount
{

    const word cellDistFieldName_; 
    const word sqrSearchDistFieldName_; 

public:

    TypeName ("markerFieldModel");

    declareRunTimeSelectionTable (
        tmp,
        markerFieldModel,
        Dictionary,
        (
            const dictionary& configDict
        ),
        (configDict)
    )

    // Constructors
    explicit markerFieldModel(const dictionary& configDict); 

    // Selectors
    static tmp<markerFieldModel> New(const dictionary& configDict);

    // Destructor
    virtual ~markerFieldModel() = default;

    // Member Functions
    virtual void calcMarkerField(volScalarField& markerField) const = 0;

    word cellDistFieldName() const
    {
        return cellDistFieldName_; 
    }

    word sqrSearchDistFieldName() const
    {
        return sqrSearchDistFieldName_; 
    }

    bool distanceWithinNarrowBand(scalar distance, scalar narrowBandWidth) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
