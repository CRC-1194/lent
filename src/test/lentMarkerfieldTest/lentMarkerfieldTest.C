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

#include "lentMarkerfieldTest.H"

namespace Foam {
namespace FrontTracking {

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
label lentMarkerfieldTest::numberInterfaceCells() const
{
    label count = 0;

    forAll(markerField_, I)
    {
        if (markerField_[I] > 0.0 && markerField_[I] < 1.0)
        {
            count++;
        }
    }

    return count;
}

void lentMarkerfieldTest::markerFieldVolumes()
{
    interfaceVolMarkerField_ = 0.0;
    frontVolMarkerField_ = 0.0;

    const fvMesh& mesh = markerField_.mesh();

    forAll(markerField_, I)
    {
        if (markerField_[I] == 0.0)
        { 
            frontVolMarkerField_ += mesh.V()[I];
        }
        else if (markerField_[I] > 0.0 && markerField_[I] < 1.0)
        {
            frontVolMarkerField_ += markerField_[I] * mesh.V()[I];
            interfaceVolMarkerField_ += markerField_[I] * mesh.V()[I];
        }
    }
}

void lentMarkerfieldTest::meshVolumes()
{
    meshVolume_ = 0.0;
    interfaceVolume_ = 0.0;

    const fvMesh& mesh = markerField_.mesh();

    forAll(markerField_, I)
    {
        if (markerField_[I] > 0.0 && markerField_[I] < 1.0)
        {
            interfaceVolume_ += mesh.V()[I];
        }

        meshVolume_ += mesh.V()[I];
    }
}

void lentMarkerfieldTest::geometricVolumes()
{
    frontVolGeometric_ = 0.0;
    interfaceVolGeometric_ = 0.0;

    // TODO: create buffer functionality in cutCellVolumeCalculator
    // to avoid duplicate, expensive computations
    scalar buffer = 0.0;

    // Iterate over all cells
    forAll(markerField_, I)
    {
        buffer = localVolCalc_.cellVolumeNegativePhase(I);

        frontVolGeometric_ += buffer;

        if (markerField_[I] > 0.0 && markerField_[I] < 1.0)
        {
            interfaceVolGeometric_ += buffer;
        }
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
lentMarkerfieldTest::lentMarkerfieldTest(const volScalarField& markerField,
                                         const triSurfaceFront& front,
                                         const dictionary& configDict)
    :
    markerField_(markerField),
    front_(front),
    configDict_(configDict),
    localVolCalc_
    (
        markerField.mesh(),
        front,
        configDict.lookup("cellDistance"),
        configDict.lookupOrDefault<word>("pointDistance","pointSignedDistance")
    )
{
    markerFieldVolumes();
    meshVolumes();
    geometricVolumes();
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

scalar lentMarkerfieldTest::interfaceVolume() const
{
    return fabs(interfaceVolMarkerField_ - interfaceVolGeometric_)
                / interfaceVolume_;
}

List<scalar> lentMarkerfieldTest::localVolume() const
{
    // Note: contrary to the other tests of this class, this 
    // function evaluates differences in terms of relative filling levels
    // rather than absolute volume differences
    List<scalar> errorSet(numberInterfaceCells(), 0.0);
    label count = 0;
    const fvMesh& mesh = markerField_.mesh();

    forAll(markerField_, I)
    {
        if (markerField_[I] > 0.0 && markerField_[I] < 1.0)
        {
            errorSet[count] = fabs(markerField_[I] -
                          localVolCalc_.cellVolumePositivePhase(I) / mesh.V()[I]);
            count++;
        }
    }

    return errorSet;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
