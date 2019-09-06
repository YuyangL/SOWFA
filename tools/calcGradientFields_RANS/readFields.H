Info<< "Reading mean velocity field U\n" << endl;
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

Info<< "Reading mean TKE field k\n" << endl;
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

Info<< "Reading mean modified pressure field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading mean temperature field T\n" << endl;
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
