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

\*---------------------------------------------------------------------------*/

#include "volPointInterpolation.H"

namespace Foam {
    namespace frontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
triSurfaceMeshDistanceCalculator::initCellSearchDistance(
    const fvMesh& mesh
)
{
    if (cellSearchDistPtr_.empty())
    {
        cellSearchDistPtr_ = autoPtr<volScalarField> (
            new volScalarField (
               IOobject (
                   "cellSearchDist", 
                   mesh.time().timeName(), 
                   mesh, 
                   IOobject::NO_READ, 
                   IOobject::NO_WRITE
               ), 
               mesh, 
               dimensionedScalar ( 
                  "zero", 
                  dimLength, 
                  0
               )
           )
        );
    }
    else
    {
        cellSearchDistPtr_->resize(mesh.cells().size()); 
    }
    cellSearchDistPtr_() = dimensionedScalar 
        (
            "GREAT", 
            dimLength, 
            GREAT
        );
}

void
triSurfaceMeshDistanceCalculator::initPointSearchDistance(
    const fvMesh& mesh
)
{
    if (pointSearchDistPtr_.empty())
    {
        pointSearchDistPtr_ = autoPtr<pointScalarField>  
        (
            new pointScalarField
            (
                IOobject
                (
                    "pointSearchDist", 
                    mesh.time().timeName(), 
                    mesh, 
                    IOobject::NO_READ, 
                    IOobject::NO_WRITE
                ),
                pointMesh(mesh), 
                0,
                "zeroGradient"
            ) 
        );
                
    }
    else
    {
        pointSearchDistPtr_->resize(mesh.points().size()); 
    }

    pointSearchDistPtr_() = GREAT;
}



void triSurfaceMeshDistanceCalculator::calcCellSearchDistance(const fvMesh& mesh)
{
    initCellSearchDistance(mesh); 

    volScalarField& cellSearchDist_ = cellSearchDistPtr_(); 

    // Sum deltaCoeffs inversed.
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs(); 

    const labelList& own = mesh.owner(); 
    const labelList& nei = mesh.neighbour(); 

    // Sum the deltaCoeffs for the internal faces.
    forAll(own, I)
    {
        cellSearchDist_[own[I]] += (1 / deltaCoeffs[I]);
        cellSearchDist_[nei[I]] += (1 / deltaCoeffs[I]);
    }

    // Sum the deltaCoeffs for the boundary faces.
    forAll(mesh.boundary(), patchI)
    {
        const fvsPatchField<scalar>& deltaCoeffsBoundary = 
            deltaCoeffs.boundaryField()[patchI];

        const labelList& faceCells =
            mesh.boundary()[patchI].faceCells();

        forAll(mesh.boundary()[patchI], faceI)
        {
            cellSearchDist_[faceCells[faceI]] += deltaCoeffsBoundary[faceI];
        }
    }

    // Correct the cell centred distance. 
    forAll(cellSearchDist_, I)
    {
        // Average the distance with the number of cell-faces. 
        cellSearchDist_[I] /=  mesh.cells()[I].size();
        // Expand the distance by the bandwidth.
        cellSearchDist_[I] *=  bandwidth_;  
    }
}

void triSurfaceMeshDistanceCalculator::calcPointSearchDistance(const fvMesh& mesh)
{
    initPointSearchDistance(mesh); 

    volPointInterpolation interpolate(mesh); 
    interpolate.interpolate(cellSearchDistPtr_(), pointSearchDistPtr_()); 
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

triSurfaceMeshDistanceCalculator::triSurfaceMeshDistanceCalculator(label bandwidth)
:
    cellsElementNearest_(),
    pointsElementNearest_(), 
    cellSearchDistPtr_(),
    pointSearchDistPtr_(),
    bandwidth_(bandwidth)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<typename Connection>
void triSurfaceMeshDistanceCalculator::calcCentresToElementsDistance(
    volScalarField& Psi, 
    const Connection& connection
) 
{
    // Get the reference to the levelSetFront.
    triSurfaceMesh front (
        IOobject
        (
        ),
        connection.front()
    );

    // Get the reference to the fvMesh
    const fvMesh& mesh = connection.mesh();

    // Compute the search distance field.
    calcCellSearchDistance(mesh); 
    
    // Get the cell centres.  
    const volVectorField& C = mesh.C(); 

    // Get the distance field reference.
    const volScalarField& cellSearchDist_ = cellSearchDistPtr_(); 

    DynamicList<pointIndexHit> cellsElementNearest_(C.size());

    mesh.time().cpuTimeIncrement();
    front.findNearest(
        C, 
        cellSearchDist_, 
        cellsElementNearest_
    );
    Info << "findNearest : " << mesh.time().cpuTimeIncrement() << endl; 

    // Create a list of the volume types: based on the cell centre, the
    // volume can be INSIDE, OUTSIDE or UNKNOWN with respect to the surface.
    List<searchableSurface::volumeType> volType;
    // Fill the list of the volume types. 
    front.getVolumeType(C, volType);

    // Get front elements.
    const List<labelledTri>& elements = front.localFaces(); 
    // Get front vertices. 
    const pointField& vertices = front.localPoints(); 

    // For all volume types. 
    forAll(volType, I)
    {
        // Get the volume type.
        searchableSurface::volumeType vT = volType[I];

        const pointIndexHit& h = cellsElementNearest_[I]; 

        if (h.hit())
        {
            // If the volume is OUTSIDE.
            if (vT == searchableSurface::OUTSIDE)
            {
                // Set the negative distance.
                Psi[I] = -Foam::mag(C[I] - h.hitPoint());
            }
            // If the volume is inside.
            else if (vT == searchableSurface::INSIDE)
            {
                // Set the positive distance.
                Psi[I] = Foam::mag(C[I] - h.hitPoint());
            }
            // FIXME 
            else // The cell is cut by the element.
            {
                // Compute the distance vector.
                //vector distance = C[I] - h.hitPoint(); 
                // Get the element.
                //const labelledTri& element = elements[h.index()];
                // Get the element normal
                //vector elementNormal = element.normal(vertices);

                // Project the distance to the element normal and set 
                // signed the distance value.
                //Psi[I] = distance & (elementNormal / mag(elementNormal));
            }

        }
    }
}

template<typename Connection>
void triSurfaceMeshDistanceCalculator::calcPointsToElementsDistance(
    scalarField& psi, 
    const Connection& connection
) 
{
    // Get the reference to the levelSetFront.
    const triSurfaceMesh& front = connection.front(); 

    // Get the reference to the fvMesh
    const fvMesh& mesh = connection.mesh();

    // Compute the search distance field.
    calcPointSearchDistance(mesh); 
    
    // Get the cell centres.  
    const pointField& points  = mesh.points(); 

    // Get the distance field reference.
    const pointScalarField& pointSearchDist_ = pointSearchDistPtr_(); 

    front.findNearest(
        points, 
        pointSearchDist_, 
        pointsElementNearest_
    );

    // Create a list of the volume types: based on the cell centre, the
    // volume can be INSIDE, OUTSIDE or UNKNOWN with respect to the surface.
    List<searchableSurface::volumeType> volType;
    // Fill the list of the volume types. 
    front.getVolumeType(points, volType);

    // Get front elements.
    const List<labelledTri>& elements = front.localFaces(); 
    // Get front vertices. 
    const pointField& vertices = front.localPoints(); 

    // For all volume types. 
    forAll(volType, I)
    {
        // Get the volume type.
        searchableSurface::volumeType vT = volType[I];

        const pointIndexHit& h = pointsElementNearest_[I]; 

        if (h.hit())
        {
            // If the volume is OUTSIDE.
            if (vT == searchableSurface::OUTSIDE)
            {
                // Set the positive distance.
                psi[I] = Foam::mag(points[I] - h.hitPoint());
            }
            // If the volume is inside.
            else if (vT == searchableSurface::INSIDE)
            {
                // Set the negative distance.
                psi[I] = -Foam::mag(points[I] - h.hitPoint());
            }
            else // The cell is cut by the element.
            {
                // Compute the distance vector.
                vector distance = points[I] - h.hitPoint(); 
                // Get the element.
                const labelledTri& element = elements[h.index()];
                // Get the element normal
                vector elementNormal = element.normal(vertices);

                // Project the distance to the element normal and set 
                // signed the distance value.
                psi[I] = distance & (elementNormal / mag(elementNormal));
            }

        }
    }


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace frontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
