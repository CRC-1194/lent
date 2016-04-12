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
    Foam::lentMarkerfieldTest

SourceFiles
    lentMarkerfieldTest.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de
    Tobias Tolle   tolle@csi.tu-darmstadt.de

Description
    This class provides all tests related to the marker field model

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

#include "fvCFD.H"

#include "lentMarkerfieldTest.H"

#include <utility>

namespace Foam
{
    namespace FrontTracking
    {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void lentMarkerfieldTest::markerfieldVolumes()
{
    phase0BulkVol_ = 0.0;
    phase0InterfaceVol_ = 0.0;
    phase1BulkVol_ = 0.0;
    phase1InterfaceVol_ = 0.0;

    const fvMesh& mesh = markerField_.mesh();

    forAll(markerField_, I)
    {
        if (markerField_[I] == 1.0)
        { 
            phase0BulkVol_ += mesh.V()[I];
        }
        else if (markerField_[I] > 0.0 && markerField_[I] < 1.0)
        {
            phase0InterfaceVol_ += markerField_[I] * mesh.V()[I];
            phase1InterfaceVol_ += (1.0 - markerField_[I]) * mesh.V()[I];
        }
        else
        {
            phase1BulkVol_ += mesh.V()[I];
        }
    }
}

void lentMarkerfieldTest::meshVolume()
{
    dimensionedScalar volume(
        "zero",
        pow(dimLength, 3),
        0.0
    );

    const fvMesh& mesh = markerField_.mesh();

    volume = sum(mesh.V());

    meshVol_ = volume.value();
}

void lentMarkerfieldTest::frontVolume()
{
    frontVolGeometric_ = 0.0;

    const List<labelledTri>& faces = front_.localFaces();
    const Field<point>& vertices = front_.localPoints();

    forAll(faces, I)
    {
        frontVolGeometric_ += tetrahedralVolume(faces[I], vertices);
    }
}

void lentMarkerfieldTest::ensureFrontPhaseIsZero()
{
    scalar phase0Vol = phase0BulkVol_ + phase0InterfaceVol_;

    if (fabs(meshVol_ - phase0Vol - frontVolGeometric_)
        < fabs(phase0Vol - frontVolGeometric_))
    {
        std::swap(phase0BulkVol_, phase1BulkVol_);
        std::swap(phase0InterfaceVol_, phase1InterfaceVol_);
    }

    frontVolMarkerField_ = phase0BulkVol_ + phase0InterfaceVol_;
}

void lentMarkerfieldTest::innerPoint()
{
    // FIXME: works only for convex shapes
    point lowerRef(GREAT, GREAT, GREAT);
    point upperRef(0.0, 0.0, 0.0);

    const Field<point>& vertices = front_.localPoints();

    forAll(vertices, I)
    {
        forAll(vertices[I], K)
        {
            if (vertices[I][K] < lowerRef[K])
            {
                lowerRef[K] = vertices[I][K];
            }

            if (vertices[I][K] > upperRef[K])
            {
                upperRef[K] = vertices[I][K];
            }
        }
    }
    
    innerPoint_ =  0.5*(lowerRef + upperRef);
}

scalar lentMarkerfieldTest::tetrahedralVolume(const labelledTri& face,
                                              const Field<point>& vertices)
{
    // use property of scalar triple product
    return fabs(1.0/6.0 * ((vertices[face[1]] - vertices[face[0]]) ^ 
                    (vertices[face[2]] - vertices[face[0]]))
                    & (innerPoint_ - vertices[face[0]]));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

lentMarkerfieldTest::lentMarkerfieldTest(const volScalarField& markerField,
                                         const triSurfaceFront& front,
                                         const dictionary& configDict)
    :
    markerField_(markerField),
    front_(front),
    configDict_(configDict)
{
    markerfieldVolumes();
    meshVolume();
    innerPoint();
    frontVolume();
    ensureFrontPhaseIsZero();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

lentMarkerfieldTest::~lentMarkerfieldTest()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool lentMarkerfieldTest::boundedness() const
{
    bool bounded = true;

    forAll(markerField_, I)
    {
        if (markerField_[I] < 0.0 || markerField_[I] > 1.0)
        {
            bounded = false;
        }
    }

    return bounded;
}

scalar lentMarkerfieldTest::globalVolume() const
{
    return fabs(frontVolMarkerField_ - frontVolGeometric_) / frontVolGeometric_;
}

scalar lentMarkerfieldTest::globalVolumeNormalized() const
{
    scalar result;

    scalar interfaceVolFraction = (phase0InterfaceVol_ + phase1InterfaceVol_)
                                    / meshVol_;
    
    result = fabs(frontVolMarkerField_ - frontVolGeometric_)
              / (frontVolGeometric_ * interfaceVolFraction);

    return result;
}

void lentMarkerfieldTest::localVolume() const
{
    cutCellVolumeCalculator localVolCalc(markerField_.mesh(), front_,
            configDict_.lookup("cellDistance"),
            configDict_.lookupOrDefault<word>("pointDistance",
                "pointSignedDistance"));

    label count = 0;
    forAll(markerField_, I)
    {
        scalar vol = localVolCalc.cellVolumeNegativePhase(I);
        if (vol < 0.0009 && vol > 0.0)
        {
            Info << "Volume = " << vol << endl;
            count++;
        }
    }
    Info << "# cells with vol fraction = " << count << endl;

    Info << "Implement me!" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
