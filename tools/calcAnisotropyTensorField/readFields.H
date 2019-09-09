Info<< "Reading Reynolds stress field uuPrime2\n" << endl;
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
