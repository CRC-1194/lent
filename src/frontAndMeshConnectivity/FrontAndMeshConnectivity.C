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
    tomislav.maric@gmx.com
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


//namespace Foam
//{
    //namespace frontTracking
    //{
        //template <class Mesh, class Front>
        //class ConnectivityListDLList; 

        //typedef FrontAndMeshConnectivity<fvMesh, levelSetFront, ConnectivityListDLList>
            //frontAndMeshConnectivity;

        //defineNamedTemplateTypeNameAndDebug(frontAndMeshConnectivity, 0); 
    //}
//}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh, class Front, 
    //template <typename Mesh, typename Front> class ConnectivityData>
    class ConnectivityData> 
Foam::frontTracking::FrontAndMeshConnectivity<Mesh, Front, ConnectivityData>::FrontAndMeshConnectivity(const Mesh& mesh, const Front& front)
    //:
        //MeshObject<Mesh, FrontAndMeshConnectivity<Mesh, Front, ConnectivityData> > (mesh), 
        //connectivity_()
{
    // //if (debug)
    //{
        Pout << "FrontAndMeshConenctivity : component ctor " << endl;
    //}
}


//Foam::frontTracking::FrontAndMeshConnectivity::FrontAndMeshConnectivity(const FrontAndMeshConnectivity&)
//:
    //baseClassName(),
    //data_()
//{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::frontTracking::FrontAndMeshConnectivity>
//Foam::frontTracking::FrontAndMeshConnectivity::New()
//{
    //return autoPtr<FrontAndMeshConnectivity>(new FrontAndMeshConnectivity);
//}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Mesh, class Front, 
    template <typename Mesh, typename Front> class ConnectivityData>
Foam::frontTracking::FrontAndMeshConnectivity<Mesh, Front, ConnectivityData>::~FrontAndMeshConnectivity()
{

    Pout << "FrontAndMeshConnectivity::dtor" << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
//
template<class Mesh, class Front, 
    template <typename Mesh, typename Front> class ConnectivityData>
const typename ConnectivityData::MapType&
Foam::frontTracking::FrontAndMeshConnectivity<Mesh, Front, ConnectivityData>::elementsToCells()
{
    return  elementsToCells_;
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//void Foam::frontTracking::FrontAndMeshConnectivity::operator=(const FrontAndMeshConnectivity& rhs)
//{
    //// Check for assignment to self
    //if (this == &rhs)
    //{
        //FatalErrorIn("Foam::frontTracking::FrontAndMeshConnectivity::operator=(const Foam::frontTracking::FrontAndMeshConnectivity&)")
            //<< "Attempted assignment to self"
            //<< abort(FatalError);
    //}
//}


// ************************************************************************* //
