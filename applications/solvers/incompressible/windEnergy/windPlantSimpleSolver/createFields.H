    // Create and read the potential temperature field
    Info << "Creating and reading potential temperature field, T..." << endl;
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
        mesh  // This makes GeometricField look in my time dir
    );

    // Create and read the modified pressure field
    Info << "Creating and reading modified pressure field, p_rgh..." << endl;
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

    // Create and read the velocity field
    Info << "Creating and reading velocity field, U..." << endl;
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

    // Create and calculate the velocity flux
    Info << "Creating and calculating velocity flux field, phi..." << endl;
    #include "createPhi.H"

    // Read the transport properties and set up the laminar (molecular) transport model
    Info << "Reading transport properties..." << endl;
    #include "readTransportProperties.H"

    // Read the atmospheric boundary layer specific properties
    Info << "Reading the atmospheric boundary layer properties..." << endl;
    #include "readABLProperties.H"

    // Read the LES anisotropy tensor field from time 0
    Info << "Reading the LES anisotropy tensor field, bij from time 0..." << endl;
    const volSymmTensorField bij
    (
        IOobject
        (
            "bij",
            runTime.timeName(0),  // Read from time 0 only
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        // If not found, it should read as 0 as below
        dimensionedSymmTensor
        (
            "bij",
            dimensionSet(0, 0, 0, 0, 0, 0, 0),  // bij is dimless
            symmTensor (0, 0, 0, 0, 0, 0)
        )
    );

    // // Create the LES and RANS blended anisotropy tensor field
    // Info << "Reading the LES and RANS blended anisotropy tensor field, bij_mix..." << endl;
    // volSymmTensorField bij_mix
    // (
    //     IOobject
    //     (
    //         "bij_mix",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     // Since NO_READ, initial values have to be provided
    //     dimensionedSymmTensor
    //     (
    //         "zero",
    //         dimensionSet(0, 0, 0, 0, 0, 0, 0),
    //         symmTensor (0, 0, 0, 0, 0, 0)
    //     )
    // );

    // Create the turbulence model (RANS, LES, or none)
    Info << "Creating turbulence model..." << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    // Create and read the wall shear stress field
    Info << "Reading the wall shear stress field, Rwall..." << endl;
    volSymmTensorField Rwall
    (
        IOobject
        (
            "Rwall",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create the Reynolds stress field
    // This is simply pre-allocating memory for this field
    Info << "Creating the Reynolds stress field R..." << endl;
    volSymmTensorField R
    (
        IOobject
        (
            "R",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->R()
        // // Since NO_READ, initial values have to be provided
        // dimensionedSymmTensor
        // (
        //     "R",
        //     Rwall.dimensions(),
        //     symmTensor (0, 0, 0, 0, 0, 0)
        // )
    );

    // // Create the LES and RANS blended Reynolds stress field
    // Info << "Creating the LES and RANS blended Reynolds stress field, Rmix..." << endl;
    // volSymmTensorField Rmix
    // (
    //     IOobject
    //     (
    //         "Rmix",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     turbulence->Rmix()
    //     // // Since NO_READ, initial values have to be provided
    //     // dimensionedSymmTensor
    //     // (
    //     //     "Rmix",
    //     //     Rwall.dimensions(),
    //     //     symmTensor (0, 0, 0, 0, 0, 0)
    //     // )
    // );

/*
    // These fields are created in kEpsilonABL.C?
    // k-epsilon turbulence model related fields
    // Create the TKE field
    Info<< "Creating the turbulent kinetic energy field, k..." << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->k()
    );

    // Create the TKE dissipation rate field
    Info<< "Creating the turbulent kinetic energy dissipation rate field, epsilon..." << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->epsilon()
    );

    // Create turbulent viscosity field
    Info<< "Creating mean SGS viscosity field, nut..." << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        turbulence->nut()
    );
*/

    // Create an object of the horizontalWindTurbineArray class if there
    // is to be a turbine array
    turbineModels::horizontalAxisWindTurbinesADM turbines(U);

    // Create Coriolis force vector
    Info << "Creating the Coriolis force vector, fCoriolis..." << endl;
    volVectorField fCoriolis
    (
        IOobject
        (
            "fCoriolis",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("fCoriolis",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
    );

    // Create and calculate the Boussinesq density field for computing buoyancy forces
    Info << "Creating kinematic (Boussinesq) density field, rhok..." << endl;
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - ( (T - TRef)/TRef )
    );

    // Create and read the turbulent thermal conductivity field
    Info << "Creating the kinematic thermal conductivity field, kappat..." << endl;
    volScalarField kappat
    (
        IOobject
        (
            "kappat",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create the wall temperature flux field
    Info << "Reading and creating the wall temperature flux field, qwall..." << endl;
    volVectorField qwall
    (
        IOobject
        (
            "qwall",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create and calculate the gravity potential field
    Info << "Creating and calculating the gravity potential field, gh and ghf..." << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    // Create and calculate the static pressure field
    Info << "Creating and calculating the static pressure field, p..." << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    // Set up the pressure reference cell information
    Info << "Setting up the pressure reference cell information..." << endl;
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("SIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
