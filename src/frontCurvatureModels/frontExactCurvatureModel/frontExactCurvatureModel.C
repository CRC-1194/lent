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
    Foam::frontExactCurvatureModel

SourceFiles
    frontExactCurvatureModel.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

    Curvature model used for validating the pressure-velocity coupling algorithm. 
    It reads an exact curvature field.

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

#include "frontExactCurvatureModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcAverage.H"
#include "fvcDiv.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(frontExactCurvatureModel, 0);
    addToRunTimeSelectionTable(frontCurvatureModel, frontExactCurvatureModel, Dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

frontExactCurvatureModel::frontExactCurvatureModel(const dictionary& configDict, const Time& runTime)
    :
        frontCurvatureModel(configDict, runTime), 
        exactCurvatureFieldName_(configDict.lookup("exactCurvatureField"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
 
frontExactCurvatureModel::~frontExactCurvatureModel() {} 

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


tmp<volScalarField> frontExactCurvatureModel::cellCurvature() const
{
    // Read the exact curvature field and forward it to the call site.
    tmp<volScalarField> exactCurvature( 
        volScalarField(
            IOobject(
                exactCurvatureFieldName_, 
                time().timeName(), 
                mesh(),
                IOobject::MUST_READ, 
                IOobject::NO_WRITE
            ),
            mesh()
        )
    );

    return exactCurvature;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
