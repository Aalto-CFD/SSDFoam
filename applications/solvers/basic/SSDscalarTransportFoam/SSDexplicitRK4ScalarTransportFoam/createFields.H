Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT
(
    transportProperties.lookup("DT")
);

#include "createPhi.H"

#include "createFvModels.H"
#include "createFvConstraints.H"



// --------------------------------- ADDED -------------------------------------

dimensionedScalar filterScope
 (
     transportProperties.lookup("filterScope")
 );

surfaceScalarField DxDyDz
(
    IOobject
    (
        "DxDyDz",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,2,0,0,0,0,0)
);

// Tfluct = (2*dx/pi)^2 * fvc::laplacian(DxDyDz,Tfluct)
// DxDyDz = SSD_gamma = (2*dx/pi)^2
DxDyDz = Foam::pow(mesh.surfaceInterpolation::deltaCoeffs() ,-2.0); // dx^2
scalar pi = constant::mathematical::pi;
// scalar filterScope = 2;
DxDyDz = (filterScope*filterScope)*DxDyDz/(4.0*pi*pi);

volScalarField Tprime
(
    IOobject
    (
        "Tprime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);
Tprime.correctBoundaryConditions();

volScalarField dT
(
    IOobject
    (
        "dT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);

volScalarField T_old
(
    IOobject
    (
        "T_old",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    T
);

volScalarField k1
(
    IOobject
    (
        "k1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,0,-1,1,0,0,0)
);

volScalarField k2
(
    IOobject
    (
        "k2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,0,-1,1,0,0,0)
);

volScalarField k3
(
    IOobject
    (
        "k3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,0,-1,1,0,0,0)
);

volScalarField k4
(
    IOobject
    (
        "k4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0,0,-1,1,0,0,0)
);
