Info<< "Reading mean velocity field UAvg\n" << endl;
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

Info<< "Reading mean resolved TKE field kResolved\n" << endl;
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

Info<< "Reading mean SGS TKE field kSGSmean\n" << endl;
volScalarField kSGSmean
(
    IOobject
    (
        "kSGSmean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading mean modified pressure field p_rghAvg\n" << endl;
volScalarField p_rghAvg
(
    IOobject
    (
        "p_rghAvg",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading mean temperature field TAvg\n" << endl;
volScalarField TAvg
(
    IOobject
    (
        "TAvg",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//     #   include "createPhi.H"

//     Info<< "Reading field T\n" << endl;
//     volScalarField T
//     (
//         IOobject
//         (
//             "T",
//             runTime.timeName(),
//             mesh,
//             IOobject::MUST_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh
//      );

//      // Read the gravitational acceleration
//      Info << "Reading the gravitational acceleration, g..." << endl;
//      uniformDimensionedVectorField g
//      (
//          IOobject
//          (
//              "g",
//              runTime.constant(),
//              mesh,
//              IOobject::MUST_READ,
//              IOobject::NO_WRITE
//          )
//      );

// Info << "Creating turbulence model..." << endl;
// singlePhaseTransportModel laminarTransport(U, phi);
// autoPtr<incompressible::turbulenceModel> turbulence
// (
//     incompressible::turbulenceModel::New(U, phi, laminarTransport)
// );

// // Turbulent Prandtl number
// dimensionedScalar Prt(laminarTransport.lookup("Prt"));
// /*
//     singlePhaseTransportModel laminarTransport(U, phi);
//     autoPtr<incompressible::RASModel> RASModel
//     (
//         incompressible::RASModel::New(U, phi, laminarTransport)
//     );
// */
