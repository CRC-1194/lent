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
#include <map>

template<class Mesh, class Front>
class ConnectivityListDLList
{
    // Private data
    List<DLList<label> > elementsToCells_; 
    List<DLList<label> > pointsToElements_;

    public: 

        ConnectivityListDLList(const Mesh& mesh, const Front& front)
            :
                elementsToCells_ (front.size()), 
                pointsToElements_(mesh.points().size())
        {

        }

        void insertElementToCell(label elementLabel, label cellLabel)
        {
            elementsToCells_[elementLabel].insert(cellLabel);
        }

        void insertPointToElement(label pointLabel, label elementLabel)
        {
            elementsToCells_[pointLabel].insert(elementLabel);
        }
};

template<class Mesh, class Front>
class ConnectivityMultiMap
{
    std::multimap<label, label> elementsToCells_;
    std::multimap<label, label> pointsToElements_;

    public:
        ConnectivityMultiMap (const Mesh& mesh, const Front& front)
            :
                elementsToCells_(),
                pointsToElements_()
        {}

        void insertElementToCell(label elementLabel, label cellLabel)
        {
            elementsToCells_.insert(std::pair<label, label> (elementLabel, cellLabel));
        }
        void insertPointToElement(label pointLabel, label elementLabel)
        {
            elementsToCells_.insert(std::pair<label, label> (pointLabel, elementLabel));
        }
};

template<class Mesh, class Front, template <typename Mesh, 
         typename Front> class ConnectivityData> 
class FrontAndMeshConnectivity
:
    public MeshObject<Mesh, 
        FrontAndMeshConnectivity<Mesh, Front, ConnectivityData> >
{
        // Private data
            //- Connectivity maps 
        ConnectivityData<Mesh, Front> connectivity_;

    public:

        FrontAndMeshConnectivity(const Mesh& mesh, const Front& front)
            :
                connectivity_(mesh, front)
        {
            Info << "FrontAndMeshConnectivity:: component ctor" << endl;

        }

};

// Encapsulate the distance between the front (vertices, elements) and the mesh
// (cells, faces, points).
template <class Mesh, class Front, template <typename Mesh, typename Front> class ConnectivityData>
class FrontAndMeshDistance
{
    public: 

        template<class MeshCellField>
        Foam::tmp<MeshCellField> cellToElementDistance(const Mesh& mesh, const Front& front)
        {
            // Get the connectivity between the front and mesh 
            FrontAndMeshConnectivity<Mesh, Front, ConnectivityData> frontAndMesh = 
                FrontAndMeshConnectivity<Mesh,Front, ConnectivityData>::New(mesh, front);
            
            // Initialize the minimal distance cell field

            // Get the elements to cells connectivity
            
            // For all elements to cells 
                // Compute the distance between each cell centre and front element
                // If the computed distance is smaller than current (minimal) 
                    // Store the computed distance to the minimal distance field
             
            // Return the distance field
        }

        template<class MeshPointField>
        Foam::tmp<MeshPointField> pointToElementDistance(const Mesh& mesh, const Front& front)
        {
            // Get the connectivity between the front and mesh 
            FrontAndMeshConnectivity<Mesh, Front, ConnectivityData> frontAndMesh = 
                FrontAndMeshConnectivity<Mesh,Front, ConnectivityData>::New(mesh, front);

            // Initialize the minimal distance point field 

            // VERSION 0
            
            // Get points to cells connectivity
            
            // Get cells to elements connectivity
            
            // For all mesh points

                // Get the cell list
                // For all cells
                    // Get the list of elements
                        // For all elements
                            // Compute the distance between the point 
                            // and the element
                            // If distance is smaller than minimal
                                // Set minimal distance

            // END VERSION 0
            
            // VERSION 1

            // Get the points to elements connectivity

            // For all points to elements
                // Get the elements of a point

                // For all elements of a point
                    // Compute the distance between the point and the element
                    // If distance is less than minimal 
                        // Update the minimal point to element distance
            
            // END VERSION 1


            // Return the pointField
            

        }


};

// Encapsulate the interpolation between the front fields and the mesh
// fields.
//template <class Mesh, class Front>
//class FrontAndMeshInterpolation
//{
    //public: 

        //// Interpolate from the front to the mesh
        //template<class FrontField, class MeshField>
        //void interpolate (MeshField& mf, const FrontField& ff)
        //{
            //// Get the connectivity between the front and mesh 
            //FrontAndMeshConnectivity<Mesh, Front> frontAndMesh = 
                //FrontAndMeshConnectivity<Mesh,Front>::New(mf.mesh(), ff.mesh());

        //}

        //// Interpoalte from the mesh to the front
        //template<class FrontField, class MeshField>
        //void interpolate (FrontField& ff, const MeshField& mf)
        //{
            //// Get the connectivity between the front and mesh 
            //FrontAndMeshConnectivity<Mesh, Front> frontAndMesh = 
                //FrontAndMeshConnectivity<Mesh,Front>::New(mf.mesh(), ff.mesh());
        //}
//};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

typedef FrontAndMeshConnectivity<fvMesh, frontTracking::levelSetFront, ConnectivityListDLList>
    frontAndMeshConnectivity;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    frontTracking::levelSetFront front;

    frontTracking::levelSetFrontVectorField vDisplacement
    (
        IOobject
        (
            "vDisplacement",
            runTime.timeName(), 
            runTime, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
        front, 
        dimensionedVector
        (
            "zero", 
            dimLength, 
            pTraits<vector>::zero
        )
    );

    volVectorField U 
    (
        IOobject
        (
            "U",
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ), 
        mesh, 
        dimensionedVector
        (
            "zero", 
            dimLength, 
            vector(1,0,0)
        )
    );
                                  

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
