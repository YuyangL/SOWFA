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
        mesh
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

    // Create the turbulence model (RANS, LES, or none)
    Info << "Creating turbulence model..." << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

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

    // Create sponge layer viscosity field
    Info << "Creating sponge layer viscosity field, spongeLayerViscosity..." << endl;
    label lengthDimension;
    if (spongeLayerType == "Rayleigh")
    {
        lengthDimension = 0;
    }
    else
    {
        lengthDimension = 2;
    }

    volScalarField spongeLayerViscosity
    (
        IOobject
        (
            "spongeLayerViscosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("spongeLayerViscosity",dimensionSet(0, lengthDimension, -1, 0, 0, 0, 0),spongeLayerViscosityTop)
    );

    forAll(mesh.cells(),cellI)
    {
        scalar z = mesh.C()[cellI][upIndex];
        spongeLayerViscosity[cellI] *= 0.5 *
            (
                1.0 - Foam::cos
                    (
                        Foam::constant::mathematical::pi *
                        max( (z - spongeLayerBaseHeight)/(spongeLayerTopHeight - spongeLayerBaseHeight) , 0.0 )
                    )
                
            );
    }

    forAll(spongeLayerViscosity.boundaryField(),i)
    {
        if ( !mesh.boundary()[i].coupled() )
        {
            forAll(spongeLayerViscosity.boundaryField()[i],j)
            {
                scalar z = mesh.boundary()[i].Cf()[j].z();
                spongeLayerViscosity.boundaryField()[i][j] *= 0.5 *
                    (
                        1.0 - Foam::cos
                            (
                                Foam::constant::mathematical::pi *
                                max( (z - spongeLayerBaseHeight)/(spongeLayerTopHeight - spongeLayerBaseHeight) , 0.0 )
                            )
                        
                    );
            }
        }
    }

    // Create sponge layer force vector
    Info << "Creating sponge layer force vector, fSponge..." << endl;
    volVectorField fSponge
    (
        IOobject
        (
            "fSponge",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("fSponge",dimensionSet(0, 1, -2, 0, 0, 0, 0),vector::zero)
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

    // Create the wall shear stress field
    Info << "Reading and creating the wall shear stress field, Rwall..." << endl;
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
      //p_rgh + rhok*gh
        p_rgh
    );


    // Set up the pressure reference cell information
    Info << "Setting up the pressure reference cell information..." << endl;
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );


    // Create and calculate the gravity potential field
    Info << "Creating and calculating the gravity potential field, gh and ghf..." << endl;
    
    vector hRef_ = vector::zero;
    if (pRefCell != -1)
    {
        hRef_ = mesh.C()[pRefCell];
    }
    reduce(hRef_,sumOp<vector>());
    dimensionedVector hRef("hRef",dimLength,hRef_);

    volScalarField gh("gh", g & (mesh.C() - hRef));
    surfaceScalarField ghf("ghf", g & (mesh.Cf() - hRef));


    // Compute the full pressure and adjust the pressure level.
    #include "adjustPressureLevel.H"
