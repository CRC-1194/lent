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
    Foam::frontSurfaceTensionForceModel

SourceFiles
    frontSurfaceTensionForceModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Interface for the front surface tension models. 

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


#ifndef frontSurfaceTensionForceModel_H
#define frontSurfaceTensionForceModel_H

#include "typeInfo.H"
#include "tmp.H"
#include "refCount.H"
#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "triSurfaceFront.H"
#include "frontCurvatureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontSurfaceTensionForceModel Declaration
\*---------------------------------------------------------------------------*/

class frontSurfaceTensionForceModel
    :
        public refCount
{

    word filterFieldName_; 

public:

    TypeName ("abstract");

    frontSurfaceTensionForceModel() = default; 
    explicit frontSurfaceTensionForceModel(const dictionary& configDict); 

    declareRunTimeSelectionTable (
        tmp,
        frontSurfaceTensionForceModel,
        Dictionary,
        (const dictionary& configDict),
        (configDict)
    );

    // Selectors
    static tmp<frontSurfaceTensionForceModel> New(const dictionary& configDict);

    // Destructor
    virtual ~frontSurfaceTensionForceModel() = default;

    const word filterFieldName() const
    {
        return filterFieldName_;  
    } 
    
    virtual tmp<surfaceScalarField> faceSurfaceTensionForce(
        const fvMesh&, 
        const triSurfaceFront&
    ) const = 0;

    virtual tmp<volVectorField> cellSurfaceTensionForce(
        const fvMesh&, 
        const triSurfaceFront&
    ) const = 0;  
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
