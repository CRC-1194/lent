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
    Foam::lentMethod

SourceFiles
    lentMethod.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description
    A facade class that simplifies the user interface for the LENT hybrid
    front-tracking level-set method by agglomerating a bunch of run-time
    selected small one-algorithm SRP classes.

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


#include "lentMethod.H"
#include "error.H"

namespace Foam  {

namespace FrontTracking {

    defineTypeNameAndDebug(lentMethod, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentMethod::lentMethod(
    const triSurfaceMesh& frontMesh,
    const fvMesh& mesh,
    word dictName
)
:
    regIOobject(
        IOobject(
            "lentMethod",
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    lentControlDict_(
        IOobject(
           dictName,
           "system",
           mesh.time(),
           IOobject::MUST_READ_IF_MODIFIED,
           IOobject::AUTO_WRITE
        )
    ),
    communicationMaps_(frontMesh, mesh),  
    frontIsReconstructed_(false),
    distanceFieldCalculatorTmp_(
       distanceFieldCalculator::New(
           lentControlDict_.subDict("distanceCalculator")
       )
    ),
    frontReconstructionModelTmp_(
       frontReconstructionModel::New(
           lentControlDict_.subDict("frontReconstructionModel")
       )
    ),
    frontReconstructorTmp_(
       frontReconstructor::New(
           lentControlDict_.subDict("frontReconstructor")
       )
    ),
    searchAlgorithmTmp_(
       frontMeshSearch::New(
           lentControlDict_.subDict("searchAlgorithm")
       )
    ),
    frontVelocityCalculatorTmp_(
       frontVelocityCalculator::New(
           lentControlDict_.subDict("frontVelocityCalculator")
       )
    ),
    frontMotionSolverTmp_(
       frontMotionSolver::New(
           lentControlDict_.subDict("frontMotionSolver")
       )
    ),
    markerFieldModelTmp_(
        markerFieldModel::New(
           lentControlDict_.subDict("markerFieldModel")
        )
    ),
    surfaceTensionForceModelTmp_(
        frontSurfaceTensionForceModel::New(
           lentControlDict_.subDict("surfaceTensionForceModel") 
        )
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void lentMethod::calcSearchDistances(
    volScalarField& searchDistanceSqr,
    pointScalarField& pointSearchDistanceSqr
)
{
    distanceFieldCalculator& distanceCalc = distanceFieldCalculatorTmp_();

    distanceCalc.calcCellSearchDistance(searchDistanceSqr);
    distanceCalc.calcPointSearchDistance(pointSearchDistanceSqr, searchDistanceSqr);
}

void lentMethod::calcSignedDistances(
    volScalarField& signedDistance,
    pointScalarField& pointSignedDistance,
    const volScalarField& searchDistanceSqr,
    const pointScalarField& pointSearchDistanceSqr,
    const triSurfaceFront& front
)
{
    // FIXME: The distance calculator wrongy re-initializes the frontMesh. TM. 
    distanceFieldCalculatorTmp_->calcCellsToFrontDistance(
        signedDistance,
        searchDistanceSqr,
        front
    );

    // FIXME: The  distance calculator wrongy re-initializes the frontMesh. TM.
    distanceFieldCalculatorTmp_->calcPointsToFrontDistance(
        pointSignedDistance,
        pointSearchDistanceSqr,
        front
    );
}

void lentMethod::calcMarkerField(volScalarField& markerField) const
{
    markerFieldModelTmp_->calcMarkerField(markerField);
}

void lentMethod::reconstructFront(
    triSurfaceFront& front,
    const volScalarField& signedDistance,
    const pointScalarField& pointSignedDistance
)
{
    if (frontReconstructionModelTmp_->reconstructionRequired(front, signedDistance))
    {
        Info << "Reconstructing front..." << endl;
        communicationMaps_.setTriangleToCell(
            frontReconstructorTmp_->reconstructFront(
                front,
                signedDistance,
                pointSignedDistance
            )
        );

        frontIsReconstructed_ = true;
        Info << "Done." << endl;
    }
}

void lentMethod::calcFrontVelocity(
    triSurfaceVectorField& frontVelocity,
    const volVectorField& U
)
{
    frontVelocityCalculatorTmp_->calcFrontVelocity(
        frontVelocity,
        U,
        communicationMaps_.triangleToCell()  // TODO: Port mesh search to the map. TM.
    );
}

void lentMethod::evolveFront(
    triSurfaceFront& front,
    const triSurfaceVectorField& frontVelocity
) const
{
    frontMotionSolverTmp_->evolveFront(
        front,
        frontVelocity
    );

    frontIsReconstructed_ = false;

    // The normals must be calculated after motion .
    calcFrontNormals(front); 
}

bool lentMethod::writeData(Ostream& os) const
{
    FatalErrorIn("lentMethod::writeData(Ostream& os)")
    << "lentMethod is not supposed to be written "
    << "regIOobject inherited to allow registry queries." << endl;

    return false;
}

void lentMethod::calcFrontNormals(triSurfaceFront& front) const
{
    auto& normals = front.storedFaceNormals(); 
    const auto& points = front.points(); 

    forAll(normals, faceI)
    {
        normals[faceI] = front[faceI].normal(points);  
        normals[faceI] /= mag(normals[faceI]) + VSMALL;
    }
}; 

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

//void lentMethod::operator=(const lentMethod& rhs)
//{
    //// Check for assignment to self
    //if (this == &rhs)
    //{
        //FatalErrorIn("lentMethod::operator=(const lentMethod&)")
            //<< "Attempted assignment to self"
            //<< abort(FatalError);
    //}

    //communicationMaps_= rhs.communicationMaps_; // FIXME Can't assign, because of regIOobject. TM.
    //lentControlDict_ = rhs.lentControlDict_;
    //distanceFieldCalculatorTmp_ = rhs.distanceFieldCalculatorTmp_;
    //frontReconstructorTmp_ = rhs.frontReconstructorTmp_;
    //markerFieldModelTmp_ = rhs.markerFieldModelTmp_;
//}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
