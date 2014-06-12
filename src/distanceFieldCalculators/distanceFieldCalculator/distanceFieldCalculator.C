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
    Foam::distanceFieldCalculator

SourceFiles
    diffuseInterfaceProperties.C

Author
    Tomislav Maric maric@csi.tu-darmstadt.de

Description

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


#include "distanceFieldCalculator.H"
#include "volPointInterpolation.H"
#include "error.H"

namespace Foam {
namespace FrontTracking {

    defineTypeNameAndDebug(distanceFieldCalculator, 0);
    defineRunTimeSelectionTable(distanceFieldCalculator, Dictionary);

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

distanceFieldCalculator::distanceFieldCalculator(
    const dictionary& configDict
)
:
    narrowBandWidth_(0)
{
    narrowBandWidth_ = readScalar(configDict.lookup("narrowBandWidth"));
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

tmp<distanceFieldCalculator>
distanceFieldCalculator::New(
   const dictionary& configDict
)
{
    const word name = configDict.lookup("type");

    // Find the constructor pointer for the model in the constructor table.
    DictionaryConstructorTable::iterator cstrIter =
        DictionaryConstructorTablePtr_->find(name);

    // If the constructor pointer is not found in the table.
    if (cstrIter == DictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn (
            "distanceFieldCalculator::New(const word& name)"
        )   << "Unknown distanceFieldCalculator type "
            << name << nl << nl
            << "Valid distanceFieldCalculators are : " << endl
            << DictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    // Construct the model and return the autoPtr to the object.
    return tmp<distanceFieldCalculator> (cstrIter()(configDict));
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
            searchDistanceSqr[faceCells[faceI]] += (1 /
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

    searchDistanceSqr.correctBoundaryConditions();
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
    pointScalarField& pointSearchDistanceSqr
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