
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

// // Read the temperature field
// Info << "Reading instantenous temperatire field, T..." << endl;
// volScalarField T
// (
//     IOobject
//     (
//         "T",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );

// // Read the gravitational acceleration
// Info << "Reading the gravitational acceleration, g..." << endl;
// uniformDimensionedVectorField g
// (
//     IOobject
//     (
//         "g",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE
//     )
// );

// // Create and calculate the velocity flux
// Info << "Creating and calculating velocity flux field, phi..." << endl;
// #include "createPhi.H"

// Info << "Creating turbulence model..." << endl;
// singlePhaseTransportModel laminarTransport(U, phi);
// autoPtr<incompressible::turbulenceModel> turbulence
// (
//     incompressible::turbulenceModel::New(U, phi, laminarTransport)
// );

/*
// If uuPrim2 exists, it's from SOWFA LES and calculate mean turbulence production
if (IOobject("uuPrime2", runTime.timeName(), mesh).headerOk())
{
    Info<< "Reading mean resolved Reynolds stress field, uuPrime2\n" << endl;
    volSymmTensorField uuPrime2
    (
        IOobject
        (
            "uuPrime2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading mean SGS viscosity field, nuSGSmean\n" << endl;
    volScalarField nuSGSmean
    (
        IOobject
        (
            "nuSGSmean",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading mean velocity field, UAvg\n" << endl;
    volVectorField UAvg
    (
        IOobject
        (
            "UAvg",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading mean resolved TKE field, kResolved\n" << endl;
    volScalarField kResolved
    (
        IOobject
        (
            "kResolved",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volSymmTensorField Rij
    (
        IOobject
        (
            "Rij",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
}
// Otherwise, it's SOWFA RANS
else
{
    // Read turbulent kinetic energy field
    Info<< "Reading (LES/ML blended) Reynolds stress field Rij..." << endl;
    volSymmTensorField Rij
    (
        IOobject
        (
            "Rij",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Read turbulent kinetic energy field
    Info<< "Reading turbulent kinetic energy field, k..." << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Read turbulent viscosity field
    Info<< "Reading turbulent viscosity field nut..." << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volSymmTensorField uuPrime2
    (
        IOobject
        (
            "uuPrime2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField nuSGSmean
    (
        IOobject
        (
            "nuSGSmean",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volVectorField UAvg
    (
        IOobject
        (
            "UAvg",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    volScalarField kResolved
    (
        IOobject
        (
            "kResolved",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
}
*/
