/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | SSDFoam: Scale-selective discretization for
   \\    /   O peration     |          OpenFOAM
    \\  /    A nd           | Website: https://github.com/moreff/SSDFoam
     \\/     M anipulation  | 2021 Ilya Morev and Ville Vuorinen
-------------------------------------------------------------------------------
License
    This file is part of SSDFoam library, derived from OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "gaussScaleSelectiveConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
const surfaceInterpolationScheme<Type>&
gaussScaleSelectiveConvectionScheme<Type>::interpScheme() const
{
    NotImplemented;
    return NullObjectRef<surfaceInterpolationScheme<Type>>();
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
gaussScaleSelectiveConvectionScheme<Type>::interpolate
(
    const surfaceScalarField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    GeometricField<Type, fvPatchField, volMesh> vfPrime =
        thiPassFilter_.ref().fvcLaplacian(DxDyDz_, vf);
    return thiResInterpScheme_().interpolate(vf+vfPrime)
         - tloResInterpScheme_().interpolate(vfPrime);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
gaussScaleSelectiveConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
gaussScaleSelectiveConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    NotImplemented;
    return NullObjectRef<tmp<fvMatrix<Type>>>();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
gaussScaleSelectiveConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
