/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Tomislav Maric, TU Darmstadt, Germany 
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

InNamespace
    Foam::fvm

Description
    Calculate the matrix for the Laplace-Beltrami of a field u on a surface
    given implicitly by a field phi.

SourceFiles
    fvmLaplacian.C

Author
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis, TU Darmstadt, Germany

\*---------------------------------------------------------------------------*/

#ifndef fvmLaplaceBeltrami_H 
#define fvmLaplaceBeltrami_H 

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "fvMatrix.H"
#include "zero.H"
#include "one.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
    template<class Type>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const zero&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const zero&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const one&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const one&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const dimensioned<GType>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const dimensioned<GType>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const GeometricField<GType, fvPatchField, volMesh>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const GeometricField<GType, fvPatchField, volMesh>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const tmp<GeometricField<GType, fvPatchField, volMesh>>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const tmp<GeometricField<GType, fvPatchField, volMesh>>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const GeometricField<GType, fvsPatchField, surfaceMesh>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const tmp<GeometricField<GType, fvsPatchField, surfaceMesh>>&,
        const GeometricField<Type, fvPatchField, volMesh>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const GeometricField<GType, fvsPatchField, surfaceMesh>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplace_beltrami
    (
        const tmp<GeometricField<GType, fvsPatchField, surfaceMesh>>&,
        const GeometricField<Type, fvPatchField, volMesh>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvmLaplaceBeltrami.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
