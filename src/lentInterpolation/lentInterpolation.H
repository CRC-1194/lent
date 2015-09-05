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
    Foam::lentInterpolation

SourceFiles
    lentInterpolation.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
        Front / Mesh communication maps.  

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


#ifndef lentInterpolation_H
#define lentInterpolation_H

#include "dictionary.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "volFields.H"
#include "lentCommunication.H"
#include "GeometricField.H"
#include "fvPatchField.H"
#include "volMesh.H"
#include "triSurfaceFrontFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class lentInterpolation Declaration
\*---------------------------------------------------------------------------*/

class lentInterpolation
{

public:

    TypeName("directTransfer");

    declareRunTimeSelectionTable (
        autoPtr,
        lentInterpolation,
        Empty,
        (
        ),
        ()
    );

    // Constructors
    lentInterpolation() = default;

    lentInterpolation(const lentInterpolation&) = default; 

    // Selectors
    static autoPtr<lentInterpolation> New(const dictionary& configDict);

    // Destructor
    virtual ~lentInterpolation() = default;

    // Member Functions
    template<typename Type>  
    void interpolate(
        DimensionedField<Type, triSurfaceFrontGeoMesh> const & frontVertexField, 
        GeometricField<Type, fvPatchField, volMesh>& vf 
    ); 

    template<typename Type>  
    void interpolate(
        DimensionedField<Type, triSurfaceFrontPointGeoMesh> const & frontField, 
        GeometricField<Type, fvPatchField, volMesh>& vf 
    ); 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lentInterpolationTemplates.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
