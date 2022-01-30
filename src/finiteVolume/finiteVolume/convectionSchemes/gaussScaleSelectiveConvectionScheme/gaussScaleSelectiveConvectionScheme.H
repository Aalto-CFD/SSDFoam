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

Class
    Foam::fv::gaussScaleSelectiveConvectionScheme

Description
    Second-order convection using face-gradients, Gauss' theorem and
    scale-selective discretization technique.

SourceFiles
    gaussScaleSelectiveConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gaussScaleSelectiveConvectionScheme_H
#define gaussScaleSelectiveConvectionScheme_H

#include "convectionScheme.H"
#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class gaussScaleSelectiveConvectionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussScaleSelectiveConvectionScheme
:
    public fv::convectionScheme<Type>
{
    // Private Data

        const fvMesh& mesh_;

        tmp<laplacianScheme<Type, scalar>> thiPassFilter_;

        scalar filterScope_;

        tmp<surfaceInterpolationScheme<Type>> tloResInterpScheme_;
        tmp<surfaceInterpolationScheme<Type>> thiResInterpScheme_;

        tmp<surfaceScalarField> DxDyDz_;



public:

    //- Runtime type information
    TypeName("GaussScaleSelective");


    // Constructors

        //- Construct from flux and interpolation scheme
        gaussScaleSelectiveConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const tmp<laplacianScheme<Type, scalar>> laplacianFiter,
            const scalar& filterScope,
            const tmp<surfaceInterpolationScheme<Type>>& loResScheme,
            const tmp<surfaceInterpolationScheme<Type>>& hiResScheme
        );

        //- Construct from flux and Istream
        gaussScaleSelectiveConvectionScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        );

        //- Disallow default bitwise copy construction
        gaussScaleSelectiveConvectionScheme(
            const gaussScaleSelectiveConvectionScheme&
        ) = delete;


    // Member Functions

        const surfaceInterpolationScheme<Type>& interpScheme() const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> interpolate
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> flux
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<fvMatrix<Type>> fvmDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        tmp<GeometricField<Type, fvPatchField, volMesh>> fvcDiv
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gaussScaleSelectiveConvectionScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gaussScaleSelectiveConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
