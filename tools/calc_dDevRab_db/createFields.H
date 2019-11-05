// Field initialization
Info << "\nCreating mean deviatoric Reynolds stress tensor field, devR..." << endl;
volSymmTensorField devR
(
    IOobject
    (
        "devR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "devR",
        dimensionSet(0, 2, -2, 0, 0, 0, 0),  // m^2/s^2
        symmTensor (123, 123, 123, 123, 123, 123)
    )
);

Info << "\nCreating mean deviatoric Reynolds stress j-gradient tensor field dDevRab_db..." << endl;
volTensorField dDevRab_db
(
    IOobject
    (
        "dDevRab_db",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor
    (
        "dDevRab_db",
        dimensionSet(0, 1, -2, 0, 0, 0, 0),  // m/s^2, a.k.a. acceleration
        tensor (123, 123, 123, 123, 123, 123, 123, 123, 123)
    )
);
