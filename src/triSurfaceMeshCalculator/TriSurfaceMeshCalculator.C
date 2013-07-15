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

//#include "volPointInterpolation.H"

#include "interpolationCellPoint.H"

namespace Foam {
    namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
TriSurfaceMeshCalculator::initCellSearchDistance(
    const fvMesh& mesh
)
{
    //if (cellSearchDistPtr_.empty())
    //{
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
               ), 
               "zeroGradient"
           )
        );
    //}
    //else
    //{
        //volScalarField cellSearchDist = cellSearchDistPtr_(); 
        //cellSearchDist.resize(mesh.cells().size()); 
        //cellSearchDist = dimensionedScalar 
        //(
            //"zero", 
            //dimLength, 
            //0 
        //);
        //cellSearchDist.boundaryField().evaluate(); 
    //}
}

void
TriSurfaceMeshCalculator::initPointSearchDistance(
    const fvMesh& mesh
)
{
    //if (pointSearchDistPtr_.empty())
    //{
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
                
    //}
    //else
    //{
        //pointSearchDistPtr_->resize(mesh.points().size()); 
    //}

    pointSearchDistPtr_() = 0;
}



void TriSurfaceMeshCalculator::calcCellSearchDistance(const fvMesh& mesh)
{
    initCellSearchDistance(mesh); 

    volScalarField& cellSearchDist_ = cellSearchDistPtr_(); 

    //Info << cellSearchDist_ << endl;

    // Sum deltaCoeffs inversed.
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs(); 

    const labelList& own = mesh.owner(); 
    const labelList& nei = mesh.neighbour(); 

    // Sum the deltaCoeffs for the internal faces.
    forAll(own, I)
    {
        cellSearchDist_[own[I]] += (1 / (deltaCoeffs[I] * deltaCoeffs[I]));
        cellSearchDist_[nei[I]] += (1 / (deltaCoeffs[I] * deltaCoeffs[I]));
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
            cellSearchDist_[faceCells[faceI]] += (1 / 
                    (deltaCoeffsBoundary[faceI] * 
                     deltaCoeffsBoundary[faceI]));
;
        }
    }

    // Correct the cell centred distance. 
    forAll(cellSearchDist_, I)
    {
        // Average the distance with the number of cell-faces. 
        cellSearchDist_[I] /=  mesh.cells()[I].size();
    }

    // Expand the distance by the bandwidth.
    cellSearchDist_ *= bandwidth_; 
    cellSearchDist_.boundaryField().evaluate(); 
}

void TriSurfaceMeshCalculator::calcPointSearchDistance(const fvMesh& mesh)
{
    initPointSearchDistance(mesh); 

    volPointInterpolation interpolate(mesh); 
    interpolate.interpolate(cellSearchDistPtr_(), pointSearchDistPtr_()); 

    // TODO: remove, debugging
    //pointSearchDistPtr_() * 2; 
}

bool TriSurfaceMeshCalculator::pointInCell
(
    const point& p, 
    label cellI, 
    const fvMesh& mesh
) const
{
    const cellList& cells = mesh.cells(); 
    const pointField& points = mesh.points(); 
    const faceList& faces = mesh.faces();
    const cell& cell = cells[cellI];
    const labelList& own = mesh.faceOwner(); 
    const surfaceVectorField& Cf = mesh.Cf(); 

    bool inside = true;

    forAll (cell, faceI)
    {
        face f = faces[cell[faceI]];

        if (! (cellI == own[cell[faceI]]))
        {
            f = f.reverseFace(); 
        }

        vector pf = p - Cf[cell[faceI]];  

        if ((pf & f.normal(points)) > 0)
        {
            inside = false;
            break;
        }
    }

    return inside;

};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

TriSurfaceMeshCalculator::TriSurfaceMeshCalculator(label bandwidth)
:
    frontMeshPtr_(), 
    cellsElementNearest_(),
    pointsElementNearest_(), 
    cellSearchDistPtr_(),
    pointSearchDistPtr_(),
    bandwidth_(bandwidth)
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <typename NarrowBandPropagation>
void TriSurfaceMeshCalculator::calcCentresToElementsDistance
(
    volScalarField& Psi, 
    const triSurfaceFront& front, 
    NarrowBandPropagation enforceNarrowBand
)
{
    Psi = dimensionedScalar
    (
        "distance", 
        dimLength, 
        GREAT
    );

    // Get the reference to the fvMesh
    const fvMesh& mesh = Psi.mesh();
    
    frontMeshPtr_ = autoPtr<triSurfaceMesh> 
    (
        new triSurfaceMesh (
            IOobject
            (
                "triSurfaceMesh", 
                "frontMesh", 
                mesh, 
                IOobject::NO_READ, 
                IOobject::NO_WRITE
            ),
            front // LSP
        )
    );

    triSurfaceMesh& frontMesh = frontMeshPtr_(); 

    // Compute the search distance field.
    Psi.time().cpuTimeIncrement(); 
    calcCellSearchDistance(mesh); 
    Info << "Compute the search distance: " 
        << Psi.time().cpuTimeIncrement() << endl; 

    // Get the cell centres.  
    const volVectorField& C = mesh.C(); 

    // Get the distance field reference.
    const volScalarField& cellSearchDist_ = cellSearchDistPtr_(); 

    mesh.time().cpuTimeIncrement();
    frontMesh.findNearest(
        C, 
        cellSearchDist_, 
        cellsElementNearest_
    );
    Info << "findNearest : " << mesh.time().cpuTimeIncrement() << endl; 

    // Create a list of the volume types: based on the cell centre, the
    List<searchableSurface::volumeType> volType;
    // Fill the list of the volume types. 
    frontMesh.getVolumeType(C, volType);

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
                //Psi[I] = -Foam::mag(C[I] - h.hitPoint());
                Psi[I] = Foam::mag(C[I] - h.hitPoint());
            }
            // If the volume is inside.
            else if (vT == searchableSurface::INSIDE)
            {
                // Set the positive distance.
                //Psi[I] = Foam::mag(C[I] - h.hitPoint());
                Psi[I] = -Foam::mag(C[I] - h.hitPoint());
            }
            // TODO: test, possibly remove 
            //else 
            //{
                //// Compute the distance vector.
                //vector distance = C[I] - h.hitPoint(); 
                //// Get the element.
                //const labelledTri& element = elements[h.index()];
                //// Get the element normal
                //vector elementNormal = element.normal(vertices);

                //// Project the distance to the element normal and set 
                //// signed the distance value.
                //Psi[I] = distance & (elementNormal / mag(elementNormal));
            //}
        }
    }

    Info << "Compute the Psi field: " 
        << Psi.time().cpuTimeIncrement() << endl;

    // Enforce the narrow band of Psi initialized with GREAT in the whole domain.
    enforceNarrowBand(Psi); 
}

//template<typename Connection, typename NarrowBandPropagation>
//void TriSurfaceMeshCalculator::calcCentresToElementsDistance
//(
    //volScalarField& Psi, 
    //Connection const & connection, 
    //NarrowBandPropagation enforceNarrowBand
//) 
//{
    //const triSurface& front = static_cast<const triSurface&> (connection.front());
    //calcCentresToElementsDistance(Psi, front, enforceNarrowBand); 
//}

template<typename Mesh, typename NarrowBandPropagation>
void TriSurfaceMeshCalculator::calcPointsToElementsDistance(
    pointScalarField& psi, 
    const triSurfaceFront& front, 
    const Mesh& mesh, 
    NarrowBandPropagation enforceNarrowBand
) 
{
    // TODO: check the initialization in createFields.
    psi = GREAT;

    // Get the reference to the fvMesh
    //const fvMesh& mesh = connection.mesh();

    // Compulistte the search distance field.
    calcPointSearchDistance(mesh); 
    
    // Get the cell centres.  
    const pointField& points  = mesh.points(); 

    // Get the distance field reference.
    const pointScalarField& pointSearchDist_ = pointSearchDistPtr_(); 

    frontMeshPtr_->findNearest(
        points, 
        pointSearchDist_, 
        pointsElementNearest_
    );

    // Create a list of the volume types: based on the cell centre, the
    // volume can be INSIDE, OUTSIDE or UNKNOWN with respect to the surface.
    List<searchableSurface::volumeType> volType;
    // Fill the list of the volume types. 
    frontMeshPtr_->getVolumeType(points, volType);

    // Get front elements.
    //const List<labelledTri>& elements = front.localFaces(); 

    // Get front vertices. 
    //const pointField& vertices = front.points(); 

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
        }
    }
    enforceNarrowBand(psi, mesh); 
}

void TriSurfaceMeshCalculator::calcFrontVelocity
(
    DynamicField<vector>& frontVelocity, 
    const triSurfaceFront& front,
    const volVectorField& U 
)
{
    // TODO: remove, optimization 
    U.time().cpuTimeIncrement(); 

    frontVelocity.resize(front.nPoints()); 
    frontVelocity = vector(0,0,0);

    interpolationCellPoint<vector> interpolation (U); 

    //const fvMesh& mesh = U.mesh(); 

    const List<labelledTri>& elements = front.localFaces(); 
    const pointField& vertices = front.points(); 
    // Update mesh cells after motion only !!! 
    const labelList& meshCells = front.meshCells(); 

    forAll (meshCells, I)
    {
        const triFace& element = elements[I]; 

        forAll (element, vertexI)
        {
            frontVelocity[element[vertexI]] = interpolation.interpolate
            (
                vertices[element[vertexI]],  
                meshCells[I] 
            );
        }
    }

    // TODO: remove, optimization 
    Info << "velocity interpolation: " << U.time().cpuTimeIncrement() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
