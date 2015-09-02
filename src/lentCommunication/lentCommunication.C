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


#include "lentCommunication.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(lentCommunication, 0);
    defineRunTimeSelectionTable(lentCommunication, FrontMesh);
    addToRunTimeSelectionTable(lentCommunication, lentCommunication, FrontMesh);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentCommunication::lentCommunication(
    const triSurfaceMesh& frontMesh, 
    const fvMesh& mesh
)
    :
        triangleToCell_()
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<lentCommunication>
lentCommunication::New(
        const dictionary& configDict, 
        const triSurfaceMesh& frontMesh, 
        const fvMesh& mesh
)
{
    const word name = configDict.lookup("type");

    FrontMeshConstructorTable::iterator cstrIter =
        FrontMeshConstructorTablePtr_->find(name);

    if (cstrIter == FrontMeshConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "lentCommunication::New(const word& name)"
        )   << "Unknown lentCommunication type "
            << name << nl << nl
            << "Valid lentCommunications are : " << endl
            << FrontMeshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return tmp<lentCommunication> (cstrIter()(frontMesh, mesh));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
