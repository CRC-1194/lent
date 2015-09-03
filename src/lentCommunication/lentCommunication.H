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
    Foam::lentCommunication

SourceFiles
    lentCommunication.C

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


#ifndef lentCommunication_H
#define lentCommunication_H

#include "dictionary.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "DynamicList.H"
#include "triSurfaceMesh.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class lentCommunication Declaration
\*---------------------------------------------------------------------------*/

class lentCommunication
    :
        public regIOobject 
{
    // Indirect triangle->cell map.
    DynamicList<label> triangleToCell_;

public:

    TypeName("listBased");

    declareRunTimeSelectionTable (
        autoPtr,
        lentCommunication,
        FrontMesh,
        (
            const triSurfaceMesh& frontMesh,  
            const fvMesh& mesh
        ),
        (frontMesh, mesh)
    );

    // Constructors
    lentCommunication(
        const triSurfaceMesh& frontMesh, 
        const fvMesh& mesh
    );

    lentCommunication(const lentCommunication&) = default; 

    // Selectors
    static autoPtr<lentCommunication> New(
        const dictionary& configDict, 
        const triSurfaceMesh& frontMesh, 
        const fvMesh& mesh
    );

    // Destructor
    virtual ~lentCommunication() = default;

    // Member Functions
    const auto& triangleToCell() const 
    {
        return triangleToCell_; 
    }

    auto& triangleToCell() // TODO: Remove after removing mesh search from the velocity interpolator. Follow through in lentMethod class. TM. 
    {
        return triangleToCell_; 
    }

    template<typename LabelList> 
    void setTriangleToCell(LabelList const& l)
    {
        triangleToCell_ = l; 
    }

    // Empty writer: regIOobject is used only for registry queries.
    bool writeData(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //