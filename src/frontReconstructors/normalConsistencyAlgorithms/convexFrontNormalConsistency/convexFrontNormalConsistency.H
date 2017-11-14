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
    Foam::convexFrontNormalConsistency

SourceFiles
    convexFrontNormalConsistency.C

Author
    Tobias Tolle    tolle@mma.tu-darmstadt.de

Description
    Ensures normal consistency by computing a geometric centre.
    ONLY applicable to CONVEX fronts.

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


#ifndef convexFrontNormalConsistency_H
#define convexFrontNormalConsistency_H

#include "normalConsistency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {
namespace FrontTracking {

/*---------------------------------------------------------------------------*\
                         Class convexFrontNormalConsistency Declaration
\*---------------------------------------------------------------------------*/

class convexFrontNormalConsistency
    :
        public normalConsistency
{

private:
    
    word orientation_;
    scalar orientationSign_;
    
    // Private Member Functions
    virtual void runNormalConsistencyAlgorithm(
        triSurfaceFront& front,
        const volScalarField& signedDistance,       // Not used by this alg.
        const pointScalarField& pointSignedDistance // Not used by this alg.
    ) const;

public:

    TypeName ("convexFrontNormal");

    // Constructors
        convexFrontNormalConsistency(const dictionary& configDict);

    // Selectors

        static tmp<convexFrontNormalConsistency> New(const dictionary& configDict);

    // Destructor
        virtual ~convexFrontNormalConsistency() {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FrontTracking

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
