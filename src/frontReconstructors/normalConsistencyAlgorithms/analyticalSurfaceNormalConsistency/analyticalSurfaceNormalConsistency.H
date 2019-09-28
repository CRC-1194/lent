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
    Foam::analyticalSurfaceNormalConsistency

SourceFiles
    analyticalSurfaceNormalConsistency.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description
    Uses an analytical surface description to ensure consistent normal
    orientation.
    Intended for testing purposes only.

    You may refer to this software as :
    //- full bibliographic data to be provided

    This code has been developed by :
        Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de (main developer)
    under the project supervision of :
        Holger Marschall <marschall@csi.tu-darmstadt.de> (group leader).
    
    Method Development and Intellectual Property :
    	Tomislav Maric maric@csi.tu-darmstadt.de, tomislav@sourceflux.de
    	Holger Marschall <marschall@csi.tu-darmstadt.de>
    	Dieter Bothe <bothe@csi.tu-darmstadt.de>

        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef analyticalSurfaceNormalConsistency_H
#define analyticalSurfaceNormalConsistency_H

#include "analyticalSurface.H"
#include "normalConsistency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class analyticalSurfaceNormalConsistency Declaration
\*---------------------------------------------------------------------------*/

class analyticalSurfaceNormalConsistency
    :
        public normalConsistency
{

private:
    
    bool reverseNormal_{false};
    tmp<analyticalSurface> surfaceTmp_;
    
    // Private Member Functions
    virtual void runNormalConsistencyAlgorithm(
        triSurfaceFront& front,
        const volScalarField& signedDistance,       // Not used by this alg.
        const pointScalarField& pointSignedDistance // Not used by this alg.
    ) const;

public:

    TypeName ("analyticalSurfaceNormal");

    // Constructors
        analyticalSurfaceNormalConsistency(const dictionary& configDict);
        analyticalSurfaceNormalConsistency(const analyticalSurface& surface, const bool reverseNormal=false);


    // Destructor
        virtual ~analyticalSurfaceNormalConsistency() {};
        
    // Member Functions
        bool makeFrontNormalsConsistent(triSurfaceFront& front) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
