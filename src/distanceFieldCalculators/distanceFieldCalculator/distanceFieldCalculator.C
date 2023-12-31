/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
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
    Foam::distanceFieldCalculator

SourceFiles
    distanceFieldCalculator.H

Authors:
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Signed-distance calculator interface.

Affiliations:
    Mathematical Modeling and Analysis Institute, Mathematics Department, 
    TU Darmstadt, Germany

Funding:
    German Research Foundation (DFG) - Project-ID 265191195 - SFB 1194

    German Research Foundation (DFG) - Project-ID MA 8465/1-1, 
    Initiation of International Collaboration 
    "Hybrid Level Set / Front Tracking methods for simulating 
    multiphase flows in geometrically complex systems"

\*---------------------------------------------------------------------------*/


#include "distanceFieldCalculator.H"
#include "volPointInterpolation.H"
#include "error.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(distanceFieldCalculator, 0);
    defineRunTimeSelectionTable(distanceFieldCalculator, Dictionary)

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distanceFieldCalculator::distanceFieldCalculator(
    const dictionary& configDict
)
:
    narrowBandWidth_(0)
{
    narrowBandWidth_ = configDict.get<label>("narrowBandWidth");
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<distanceFieldCalculator>
distanceFieldCalculator::New(
   const dictionary& configDict
)
{
    const word name = configDict.get<word>("type");

    // Find the constructor pointer for the model in the constructor table.
    auto* ctorPtr = DictionaryConstructorTable(name);

    // If the constructor pointer is not found in the table.
    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            configDict,
            "distanceFieldCalculator",
            name,
            *DictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    // Construct the model and return the autoPtr to the object.
    return tmp<distanceFieldCalculator> (ctorPtr(configDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

distanceFieldCalculator::~distanceFieldCalculator()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void distanceFieldCalculator::calcCellSearchDistance(
    volScalarField& searchDistanceSqr
)
{
    const fvMesh& mesh = searchDistanceSqr.mesh();

    // Sum deltaCoeffs inversed.
    const surfaceScalarField& deltaCoeffs = mesh.deltaCoeffs();

    const labelList& own = mesh.owner();
    const labelList& nei = mesh.neighbour();

    // Sum the deltaCoeffs for the internal faces.
    forAll(own, I)
    {
        searchDistanceSqr[own[I]] += (1 / (deltaCoeffs[I] * deltaCoeffs[I]));
        searchDistanceSqr[nei[I]] += (1 / (deltaCoeffs[I] * deltaCoeffs[I]));
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
            searchDistanceSqr[faceCells[faceI]] += (2 /
                    (deltaCoeffsBoundary[faceI] *
                     deltaCoeffsBoundary[faceI]));
        }
    }

    // Correct the cell centered distance.
    forAll(searchDistanceSqr, I)
    {
        // Average the distance with the number of cell-faces.
        searchDistanceSqr[I] /=  mesh.cells()[I].size();
    }

    // Expand the distance by the bandwidth.
    searchDistanceSqr == searchDistanceSqr * narrowBandWidth_ * narrowBandWidth_;
}

void distanceFieldCalculator::calcPointSearchDistance(
    pointScalarField& pointSearchDistanceSqr,
    const volScalarField& searchDistanceSqr
)
{
    const fvMesh& mesh = searchDistanceSqr.mesh();

    volPointInterpolation ip(mesh);

    pointSearchDistanceSqr.resize(mesh.nPoints());

    //pointSearchDistanceSqr = ip.interpolate(searchDistanceSqr, "fixedValue", true);
    ip.interpolate(searchDistanceSqr, pointSearchDistanceSqr);
}

void distanceFieldCalculator::calcPointSearchDistance(
    pointScalarField&
)
{
    notImplemented("distanceFieldCalculator::calcPointSearchDistance(pointScalarField&");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
