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
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeFvConvectionScheme(gaussScaleSelectiveConvectionScheme)

// ************************************************************************* //
