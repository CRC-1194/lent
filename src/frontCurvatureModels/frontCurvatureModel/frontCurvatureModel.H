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
    Foam::frontCurvatureModel

SourceFiles
    frontCurvatureModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    Interface for the front curvature models. 

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


#ifndef frontCurvatureModel_H
#define frontCurvatureModel_H

#include "typeInfo.H"
#include "tmp.H"
#include "refCount.H"
#include "volFieldsFwd.H"
#include "triSurfaceFront.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontCurvatureModel Declaration
\*---------------------------------------------------------------------------*/

class frontCurvatureModel
    :
        public refCount
{

    const Time& runTime_; 

    const word meshName_; 
    const fvMesh& mesh_; 

    const word inputFieldName_; 
    const volScalarField& inputField_; 

    // Stabilizator value for the gradient to avoid division with zero.
    const scalar delta_; 

    frontCurvatureModel(const frontCurvatureModel&); 
    void operator=(const frontCurvatureModel&);

public:

    TypeName ("base");

    declareRunTimeSelectionTable (
        tmp,
        frontCurvatureModel,
        Dictionary,
        (
            const dictionary& configDict, 
            const Time& runTime
        ),
        (configDict, runTime)
    );

    // Constructors

        frontCurvatureModel(const dictionary& configDict, const Time& runTime);

    // Selectors

        static tmp<frontCurvatureModel> New(const dictionary& configDict, const Time& runTime);

    // Destructor
        virtual ~frontCurvatureModel();

    // Member Functions
        virtual tmp<volScalarField> cellCurvature() const; 

        //virtual tmp<surfaceScalarField> faceCurvature() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
