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

Application
    protFrontAndMeshConnectivity

Description
    Prototyping (test driven development) for the connectivity between the 
    front and the mesh.

Author
    Tomislav Maric
    maric@csi.tu-darmstadt.de
    Mathematical Modelling and Analysis Group 
    Center of Smart Interfaces
    TU Darmstadt
    Germany


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "levelSetFront.H"
#include "levelSetFrontFields.H"
#include "MeshObject.H"
#include "DLList.H"
#include <list>


class MapMultiMap
{
};

template<class Mesh, class Front, class ConnectivityData = MapMultiMap>
class FrontAndMeshConnectivity
:
    public MeshObject<Mesh, FrontAndMeshConnectivity<Mesh, Front> >
{
    // Private data
    ConnectivityData data_;
        
    public:

        FrontAndMeshConnectivity(const Mesh& mesh, const Front& front)
        {
            Info << "FrontAndMeshConnectivity:: component ctor" << endl;
        }

        // Update the connectivity information using the new front topology.
        void update();
};

// Encapsulate the distance between the front (vertices, elements) and the mesh
// (cells, faces, points).
template <class Mesh, class Front>
class FrontAndMeshDistance
{
    // Private data
        //- Cell to element distance field
        
    public: 

        template<typename FrontField, typename MeshField>
        Foam::tmp<MeshField> cellElementDistance(const Mesh& mesh, const Front& front)
        {
            // Get the connectivity between the front and mesh 
            FrontAndMeshConnectivity<Mesh, Front> frontAndMesh = 
                FrontAndMeshConnectivity<Mesh,Front>::New(mesh, front);
        }
};

// Encapsulate the interpolation between the front fields and the mesh
// fields.
template <class Mesh, class Front>
class FrontAndMeshInterpolation
{
    public: 

        template<class FrontField, class MeshField>
        void interpolate (MeshField& mf, const FrontField& ff)
        {
            // Get the connectivity between the front and mesh 
            FrontAndMeshConnectivity<Mesh, Front> frontAndMesh = 
                FrontAndMeshConnectivity<Mesh,Front>::New(mf.mesh(), ff.mesh());
        }

        template<class FrontField, class MeshField>
        void interpolate (FrontField& ff, const MeshField& mf)
        {
            // Get the connectivity between the front and mesh 
            FrontAndMeshConnectivity<Mesh, Front> frontAndMesh = 
                FrontAndMeshConnectivity<Mesh,Front>::New(mf.mesh(), ff.mesh());
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    frontTracking::levelSetFront front;

    frontTracking::levelSetFrontScalarField vCurvature
    (
        IOobject
        (
            "vCurvature",
            runTime.timeName(), 
            runTime, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
        front, 
        dimensionedScalar
        (
            "zero", 
            dimLength, 
            pTraits<scalar>::zero
        )
    );

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
