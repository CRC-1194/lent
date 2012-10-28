/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany


\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{
        defineNamedTemplateTypeNameAndDebug(frontTracking::frontAndMeshConnectivity<fvMesh>, 0); 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh>
Foam::frontTracking::frontAndMeshConnectivity<Mesh>::frontAndMeshConnectivity(const Mesh& mesh, const levelSetFront& front)
    :
        MeshObject<Mesh, frontAndMeshConnectivity<Mesh> > (mesh)
{
    Info << "frontAndMeshConnectivity:: component ctor " << endl;
}


//Foam::frontTracking::frontAndMeshConnectivity::frontAndMeshConnectivity(const frontAndMeshConnectivity&)
//:
    //baseClassName(),
    //data_()
//{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::frontTracking::frontAndMeshConnectivity>
//Foam::frontTracking::frontAndMeshConnectivity::New()
//{
    //return autoPtr<frontAndMeshConnectivity>(new frontAndMeshConnectivity);
//}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Mesh>
Foam::frontTracking::frontAndMeshConnectivity<Mesh>::~frontAndMeshConnectivity()
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//void Foam::frontTracking::frontAndMeshConnectivity::operator=(const frontAndMeshConnectivity& rhs)
//{
    //// Check for assignment to self
    //if (this == &rhs)
    //{
        //FatalErrorIn("Foam::frontTracking::frontAndMeshConnectivity::operator=(const Foam::frontTracking::frontAndMeshConnectivity&)")
            //<< "Attempted assignment to self"
            //<< abort(FatalError);
    //}
//}


// ************************************************************************* //
