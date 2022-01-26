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

Author
    Tomislav Maric, maric@mma.tu-darmstadt.de
    Mathematical Modeling and Analysis, TU Darmstadt, Germany

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrix.H"
#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return fvm::laplace_beltrami(Gamma, vf, name);
}


template<class Type>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "1",
            vf.time().constant(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return fvm::laplace_beltrami
    (
        Gamma,
        vf,
        "laplace_beltrami(" + vf.name() + ')'
    );
}


template<class Type>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const zero&,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const zero&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<fvMatrix<Type>>
    (
        new fvMatrix<Type>(vf, dimensionSet(0, 0, -2, 0, 0))
    );
}


template<class Type>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const one&,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fvm::laplace_beltrami(vf, name);
}


template<class Type>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const one&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvm::laplace_beltrami(vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    const GeometricField<GType, fvsPatchField, surfaceMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fvm::laplace_beltrami(Gamma, vf, name);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const GeometricField<GType, fvsPatchField, surfaceMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fvm::laplace_beltrami(Gamma, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    ).ref().fvmLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const tmp<GeometricField<GType, fvPatchField, volMesh>>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> Laplacian(fvm::laplace_beltrami(tgamma(), vf, name));
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvm::laplace_beltrami
    (
        gamma,
        vf,
        "laplace_beltrami(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const tmp<GeometricField<GType, fvPatchField, volMesh>>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> Laplacian(fvm::laplace_beltrami(tgamma(), vf));
    tgamma.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    ).ref().fvmLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh>>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> tLaplacian = fvm::laplace_beltrami(tgamma(), vf, name);
    tgamma.clear();
    return tLaplacian;
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvm::laplace_beltrami
    (
        gamma,
        vf,
        "laplace_beltrami(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<fvMatrix<Type>>
laplace_beltrami
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh>>& tGamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm(fvm::laplace_beltrami(tGamma(), vf));
    tGamma.clear();
    return tfvm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
