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

Application
    explicitRK4scalarTransportFoam

Description
    Solves the steady or transient transport equation for a passive scalar using
    4th order Runge-Kutta explicit time-stepping.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        // RK4: step 1
        T_old = T;
        k1 = (
           - fvc::div(phi, T)
           + fvc::laplacian(DT, T)
        );
        T = T_old + 0.5*k1*runTime.deltaT();
        // RK4: step 2
        k2 = (
           - fvc::div(phi, T)
           + fvc::laplacian(DT, T)
        );
        T = T_old + 0.5*k2*runTime.deltaT();
        // RK3: step 3
        k3 = (
           - fvc::div(phi, T)
           + fvc::laplacian(DT, T)
        );
        T = T_old + k3*runTime.deltaT();
        // RK4: step 4
        k4 = (
           - fvc::div(phi, T)
           + fvc::laplacian(DT, T)
        );
        // RK4: result
        T = T_old + (1./6.)*runTime.deltaT()*(k1 + 2*k2 + 2*k3 + k4);


        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
