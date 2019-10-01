// Read the velocity field
Info << "Reading instantenous velocity field, U..." << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Read the temperature field
Info << "Reading instantenous temperatire field, T..." << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

// Read the gravitational acceleration
Info << "Reading the gravitational acceleration, g..." << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Create and calculate the velocity flux
Info << "Creating and calculating velocity flux field, phi..." << endl;
#   include "createPhi.H"

Info << "Creating turbulence model..." << endl;
singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// // Read mean SGS turbulent viscosity field
// Info<< "Reading mean SGS viscosity field, nuSGSmean..." << endl;
// volScalarField nuSGSmean
// (
//     IOobject
//     (
//         "nuSGSmean",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );

// Read moleculer viscosity field
Info<< "Reading moleculer viscosity field, nu..." << endl;
volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    turbulence->nu()
);
