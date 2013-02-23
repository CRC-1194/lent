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

\*---------------------------------------------------------------------------*/

#include "frontTrackingCalculator.H"

namespace Foam { 
    namespace frontTracking {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//void frontTrackingCalculator::resetDistanceField(volScalarField& psi)
//{
    //// Reset the field to GREAT. 
    //psi = dimensionedScalar ("signedDistance", dimLength, GREAT); 
//}

//scalar frontTrackingCalculator::pointToElementDistance ( 
     //const point& inputPoint, 
     //label elementI
//)
//{
    //// Get the front. 
    //const levelSetFront& front = connectivity_.front(); 
   
    //// Get the front data.
    //const pointField& frontPoints = front.points(); 
    //const vectorField& elementNormals = front.faceNormals(); 
    //const List<labelledTri>& frontFaces = front.localFaces(); 
    //const labelledTri element = frontFaces[elementI];
    
    //// Compute the distance between the point and the element. 
    //scalar pointToElementDist = GREAT; 

    //// Project the point onto the element plane. 
    //plane elementPlane (frontPoints[element[0]], elementNormals[elementI]);

    //// Project the point to the plane.  
    //point inputPointProjected = elementPlane.nearestPoint(inputPoint); 

    //// If the point lies in the element
        //// Compute the signed distance
    //// Else 
        //// Get element edges
        //// Initialize the minimal absolute distance
        //// Initialize the signed distance  
        //// For all element edges
            //// Compute the distance to edge
            //// Compute the squared distance
            //// If squared distance to edge is smaller than minimal 
                //// Update the minimal distance
                //// Update the signed distance

        //// Return the signed distance corresponding to the minimal squared distance
        
    //return pointToElementDist; 
//}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class MeshAndFrontConnection>
frontTrackingCalculator<MeshAndFrontConnection>::frontTrackingCalculator (
    //const fvMesh& mesh, 
    //const levelSetFront& front
)
//: 
    // Initialize the calculators
    //connectivity_(mesh, front)
{
    Pout << " frontTrackingCalculator::frontTrackingCalculator ( "
        << "\n    const fvMesh& mesh, " 
        << "\n    const levelSetFront& front " << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

//frontTrackingCalculator::~frontTrackingCalculator()
//{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class MeshAndFrontConnection>
void frontTrackingCalculator<MeshAndFrontConnection>::calcCentresToElementsDistance (
    volScalarField& Psi, 
    const MeshAndFrontConnection& connectivity
)
{
    // Delegate distance claculation to the RTS distance claculator.
    distanceCalculator_.calcCentresToElementsDistance(Psi, connectivity); 

    // Reset the distance before calculation to GREAT

    // For all cells-elements
        // For all the cell-elements
            // Compute the distance between the cell and the element


    //Pout << "frontTrackingCalculator::calcDistanceField(volScalarField& psi)" 
        //<< endl;

    //// Re-set the field to infinite distance
    //resetDistanceField(psi); 

    //// Get the nearest elements to cells. 
    //const List<pointIndexHit> & nearest = connectivity_.cellsToElementsNearest(); 

    //// Get the cell centres.  
    //const volVectorField& C = connectivity_.mesh().C(); 
    
    //// Get the reference to the front.
    //const levelSetFront& front = connectivity_.front(); 

    //// Create a list of the volume types: based on the cell centre, the
    //// volume can be INSIDE, OUTSIDE or UNKNOWN with respect to the surface.
    //List<searchableSurface::volumeType> volType;
    //// Fill the list of the volume types. 
    //front.getVolumeType(C, volType);

    //// Get front elements.
    //const List<labelledTri>& elements = front.localFaces(); 
    //// Get front vertices. 
    //const pointField& vertices = front.localPoints(); 

    //// For all volume types. 
    //forAll(volType, I)
    //{
        //// Get the volume type.
        //searchableSurface::volumeType vT = volType[I];

        //const pointIndexHit& h = nearest[I]; 

        //if (h.hit())
        //{
            //// If the volume is OUTSIDE.
            //if (vT == searchableSurface::OUTSIDE)
            //{
                //// Set the positive distance.
                //psi[I] = Foam::mag(C[I] - h.hitPoint());
            //}
            //// If the volume is inside.
            //else if (vT == searchableSurface::INSIDE)
            //{
                //// Set the negative distance.
                //psi[I] = -Foam::mag(C[I] - h.hitPoint());
            //}
            //else // The cell is cut by the element.
            //{
                //// Compute the distance vector.
                //vector distance = C[I] - h.hitPoint(); 
                //// Get the element.
                //const labelledTri& element = elements[h.index()];
                //// Get the element normal
                //vector elementNormal = element.normal(vertices);

                //// Project the distance to the element normal and set 
                //// signed the distance value.
                //psi[I] = distance & (elementNormal / mag(elementNormal));
            //}

        //}
    //}
}

template<class MeshAndFrontConnection>
void frontTrackingCalculator<MeshAndFrontConnection>::calcPointsToElementsDistance (
    scalarField& psi, 
    const MeshAndFrontConnection& connectivity
)
{
    distanceCalculator_.calcPointsToElementsDistance(psi, connectivity);
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//void frontTrackingCalculator::operator=(const frontTrackingCalculator& rhs)
//{
    //// Check for assignment to self
    //if (this == &rhs)
    //{
        //FatalErrorIn("frontTrackingCalculator::operator=(const frontTrackingCalculator&)")
            //<< "Attempted assignment to self"
            //<< abort(FatalError);
    //}
//}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
