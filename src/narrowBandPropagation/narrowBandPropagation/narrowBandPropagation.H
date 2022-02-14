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
    Foam::narrowBandPropagation

SourceFiles
    narrowBandPropagation.C

Authors
    Tomislav Maric (maric@mma.tu-darmstadt.de)

Description
    Class for propagating the switch in the sign of the sign distance field.
    Obsolete (SMCI/SMCA): 
        Tolle, Tobias, Dirk Gründing, Dieter Bothe, and Tomislav Marić.
        "triSurfaceImmersion: Computing volume fractions and signed distances
        from triangulated surfaces immersed in unstructured meshes." Computer
        Physics Communications 273 (2022): 108249.


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


#ifndef narrowBandPropagation_H
#define narrowBandPropagation_H

#include "refCount.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FrontTracking
{

/*---------------------------------------------------------------------------*\
                         Class narrowBandPropagation Declaration
\*---------------------------------------------------------------------------*/

class narrowBandPropagation
:
    public refCount
{

public:

    TypeName("none");

    declareRunTimeSelectionTable (
        tmp,
        narrowBandPropagation,
        Dictionary,
        (
            const dictionary& configDict
        ),
        (configDict)
    )

    // Constructors

        narrowBandPropagation(const dictionary& configDict);

    // Destructor
        virtual ~narrowBandPropagation();

    // Selectors

        static tmp<narrowBandPropagation> New(const dictionary& configDict);

    // Member Functions

        virtual void ensureNarrowBand(
            volScalarField&,
            scalar
        ) const 
        {} 

        virtual void ensureNarrowBand(
            pointScalarField&,
            scalar
        ) const
        {} 
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
