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
    Foam::frontCurvatureMeyer

SourceFiles
    frontCurvatureMeyer.C

Author
    Tobias Tolle bt@lefou-familie.org

Contributors 
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

    Curvature model that averages the input field and computes the curvature
    using the divergence of the normal field.

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


#ifndef frontCurvatureMeyer_H
#define frontCurvatureMeyer_H

#include "frontCurvatureModel.H"
#include "lentInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class frontCurvatureMeyer Declaration
\*---------------------------------------------------------------------------*/

class frontCurvatureMeyer
    :
        public frontCurvatureModel 
{

private:
    
    mutable tmp<triSurfaceFrontPointVectorField> curvatureNormalTmp_;

    void initializeCurvatureNormal(const fvMesh&, const triSurfaceFront&) const;


public:

    TypeName ("Meyer");

    // Constructors
    frontCurvatureMeyer(const dictionary& configDict);

    virtual ~frontCurvatureMeyer() = default;  

    // Member functions
    labelList orderVertices(labelledTri& tri, label V) const; 

    scalar getAngle(vector& a, vector& b) const;

    scalar cot(vector a, vector b) const;

    scalar areaHeron(vector& A, vector& B, vector& C) const;

    virtual tmp<volScalarField> cellCurvature(const fvMesh&, const triSurfaceFront&) const; 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
