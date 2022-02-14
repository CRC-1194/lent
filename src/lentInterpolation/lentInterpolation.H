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
    Foam::lentInterpolation

SourceFiles
    lentInterpolation.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    LENT interpolation interface
        
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


#ifndef lentInterpolation_H
#define lentInterpolation_H

#include "dictionary.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "volFields.H"
#include "lentCommunication.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "triSurfaceFrontFields.H" 
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class lentInterpolation Declaration
\*---------------------------------------------------------------------------*/

class lentInterpolation
{
    word cellToVertexInterpType_; 

public:

    TypeName("directTransfer");

    declareRunTimeSelectionTable (
        autoPtr,
        lentInterpolation,
        Empty,
        (
        ),
        ()
    )

    // Constructors
    lentInterpolation(); 

    lentInterpolation(const lentInterpolation&) = default; 

    lentInterpolation(const dictionary& dict); 

    // Selectors
    static autoPtr<lentInterpolation> New(const dictionary& configDict);

    // Destructor
    virtual ~lentInterpolation() = default;

    // Member Functions
    
    // Interpolate front vertex field onto cell centered volume field.  
    template<typename Type>  
    void interpolate(
        const DimensionedField<Type, triSurfaceFrontPointGeoMesh>& vertexField, 
        GeometricField<Type, fvPatchField, volMesh>& cellField 
    );  

    // Interpolate tmp front vertex field onto cell centered volume field.  
    template<typename Type>  
    void interpolate(
        tmp<DimensionedField<Type, triSurfaceFrontPointGeoMesh>> vertexFieldTmp, 
        GeometricField<Type, fvPatchField, volMesh>& cellField 
    );  

    // Interpolate cell centered volume field onto front vertex field.  
    template<typename Type>  
    void interpolate(
        const GeometricField<Type, fvPatchField, volMesh>& cellField, 
        DimensionedField<Type, triSurfaceFrontPointGeoMesh>& vertexField 
    );  

    // Interpolate tmp cell centered volume field onto front vertex field.  
    template<typename Type>  
    void interpolate(
        tmp<GeometricField<Type, fvPatchField, volMesh>> cellFieldTmp, 
        DimensionedField<Type, triSurfaceFrontPointGeoMesh>& vertexField
    );  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
    #include "lentInterpolationTemplates.H" 
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
