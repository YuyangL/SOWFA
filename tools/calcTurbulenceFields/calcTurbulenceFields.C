/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    calcTurbulenceFields

Description
    Calculate and create the SGS dissipation rate epsilon field as well resolved and total epsilon field.
    Furthermore, if U' exists calculated resolved and total turbulent kinetic energy field kResolved and kTotal

Source files:
    readFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info << "Time = " << runTime.timeName() << endl;

        #include "readFields.H"

        // If epsilonResolved is not calculated already, calculate it
        if (!IOobject("epsilonResolved", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nRetrieving mean resolved dissipation rate field, epsilonResolved..." << endl;
            volScalarField epsilonResolved
            (
                IOobject
                (
                    "epsilonResolved",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                2.*nu*(symm(fvc::grad(UAvg)) && symm(fvc::grad(UAvg)))
            );

            Info << "\nWriting mean resolved dissipation rate field, epsilonResolved..." << endl;
            epsilonResolved.write();

            scalar epsresolved_min = min(epsilonResolved).value();
            scalar epsresolved_max = max(epsilonResolved).value();
            Info << "Min epsilonResolved is " << epsresolved_min << "; max is " << epsresolved_max << endl;
        }
        else
        {
            Info << "\nMean resolved dissipation rate field epsilonResolved already exists!" << endl;
        }

        // If epsilonTotal is not calculated already, calculate it
        if (!IOobject("epsilonTotal", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nReading mean SGS dissipation rate field epsilonSGSmean..." << endl;
            volScalarField epsilonSGSmean
            (
                IOobject
                (
                    "epsilonSGSmean",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            Info << "\nRetrieving mean total dissipation rate field, epsilonTotal..." << endl;
            volScalarField epsilonTotal
            (
                IOobject
                (
                    "epsilonTotal",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                2.*nu*(symm(fvc::grad(UAvg)) && symm(fvc::grad(UAvg))) + epsilonSGSmean
            );

            Info << "\nWriting mean total dissipation rate field, epsilonTotal..." << endl;
            epsilonTotal.write();

            scalar epstot_min = min(epsilonTotal).value();
            scalar epstot_max = max(epsilonTotal).value();
            Info << "Min epsilonTotal is " << epstot_min << "; max is " << epstot_max << endl;
        }
        else
        {
            Info << "\nMean total dissipation rate field epsilonTotal already exists!" << endl;
        }

        // If U' exists, calculate resolved and total TKE as well
        if (IOobject("uuPrime2", runTime.timeName(), mesh).headerOk())
        {
            Info<< "Reading mean resolved turbulent stress field, uuPrime2\n" << endl;
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

            // If kResolved is not calculated already, calculate it
            if (!IOobject("kResolved", runTime.timeName(), mesh).headerOk())
            {
                Info << "\nRetrieving mean resolved turbulent kinetic energy field, kResolved..." << endl;
                volScalarField kResolved
                (
                    IOobject
                    (
                        "kResolved",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    0.5*(tr(uuPrime2))
                );

                Info << "\nWriting mean resolved turbulent kinetic energy field, kResolved..." << endl;
                kResolved.write();

                scalar kresolved_min = min(kResolved).value();
                scalar kresolved_max = max(kResolved).value();
                Info << "Min kResolved is " << kresolved_min << "; max is " << kresolved_max << endl;
            }
            else
            {
                Info << "\nMean resolved turbulent kinetic energy field kResolved already exists!" << endl;
            }

            // If kTotal is not calculated already, calculate it
            if (!IOobject("kTotal", runTime.timeName(), mesh).headerOk())
            {
                Info<< "Reading mean SGS turbulent kinetic energy field, kSGSmean\n" << endl;
                volScalarField kSGSmean
                (
                    IOobject
                    (
                        "kSGSmean",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh
                );

                Info << "\nRetrieving total turbulent kinetic energy field, kTotal..." << endl;
                volScalarField kTotal
                (
                    IOobject
                    (
                        "kTotal",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    0.5*(tr(uuPrime2)) + kSGSmean
                );

                Info << "\nWriting mean total turbulent kinetic energy field, kTotal..." << endl;
                kTotal.write();

                scalar ktot_min = min(kTotal).value();
                scalar ktot_max = max(kTotal).value();
                Info << "Min kTotal is " << ktot_min << "; max is " << ktot_max << endl;
            }
            else
            {
                Info << "\nMean total turbulent kinetic energy field kTotal already exists!" << endl;
            }
        }
        else
        {
            Info << "\nMean resolved turbulent stress field uuPrime2 not found! Not calculating mean resolved and total turbulent kinetic energy kResolved and kTotal!" << endl;
        }

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
