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

        if (!IOobject("epsilonSGS", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nRetrieving SGS dissipation rate field, epsilonSGS..." << endl;
            volScalarField epsilonSGS
            (
                IOobject
                (
                    "epsilonSGS",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                turbulence->epsilon()
            );

            Info << "\nWriting SGS dissipation rate field, epsilonSGS..." << endl;
            epsilonSGS.write();

            scalar epssgs_min = min(epsilonSGS).value();
            scalar epssgs_max = max(epsilonSGS).value();
            Info << "Min epsilonSGS is " << epssgs_min << "; max is " << epssgs_max << endl;
        }
        else
        {
            Info << "\nSGS dissipation rate field epsilonSGS already exists!" << endl;
        }

        // If epsilonResolved is not calculated already, calculate it
        if (!IOobject("epsilonResolved", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nRetrieving resolved dissipation rate field, epsilonResolved..." << endl;
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
                2.*nu*(symm(fvc::grad(U)) && symm(fvc::grad(U)))
            );

            Info << "\nWriting resolved dissipation rate field, epsilonResolved..." << endl;
            epsilonResolved.write();

            scalar epsresolved_min = min(epsilonResolved).value();
            scalar epsresolved_max = max(epsilonResolved).value();
            Info << "Min epsilonResolved is " << epsresolved_min << "; max is " << epsresolved_max << endl;
        }
        else
        {
            Info << "\nResolved dissipation rate field epsilonResolved already exists!" << endl;
        }

        // If epsilonTotal is not calculated already, calculate it
        if (!IOobject("epsilonTotal", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nRetrieving total dissipation rate field, epsilonTotal..." << endl;
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
                2.*nu*(symm(fvc::grad(U)) && symm(fvc::grad(U))) + turbulence->epsilon()
            );

            Info << "\nWriting total dissipation rate field, epsilonTotal..." << endl;
            epsilonTotal.write();

            scalar epstot_min = min(epsilonTotal).value();
            scalar epstot_max = max(epsilonTotal).value();
            Info << "Min epsTotal is " << epstot_min << "; max is " << epstot_max << endl;
        }
        else
        {
            Info << "\nTotal dissipation rate field epsilonTotal already exists!" << endl;
        }

        // If U' exists, calculate resolved and total TKE as well
        if (IOobject("Uprime", runTime.timeName(), mesh).headerOk())
        {
            Info<< "Reading fluctuating velocity field, Uprime\n" << endl;
            volVectorField Uprime
            (
                IOobject
                (
                    "Uprime",
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
                Info << "\nRetrieving resolved turbulent kinetic energy field, kResolved..." << endl;
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
                    0.5*(tr(sqr(Uprime)))
                );

                Info << "\nWriting resolved turbulent kinetic energy field, kResolved..." << endl;
                kResolved.write();

                scalar kresolved_min = min(kResolved).value();
                scalar kresolved_max = max(kResolved).value();
                Info << "Min kResolved is " << kresolved_min << "; max is " << kresolved_max << endl;
            }
            else
            {
                Info << "\nResolved turbulent kinetic energy field, kResolved already exists!" << endl;
            }

            // If kTotal is not calculated already, calculate it
            if (!IOobject("kTotal", runTime.timeName(), mesh).headerOk())
            {
                Info<< "Reading SGS turbulent kinetic energy field, k\n" << endl;
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
                    0.5*(tr(sqr(Uprime))) + k
                );

                Info << "\nWriting total turbulent kinetic energy field, kTotal..." << endl;
                kTotal.write();

                scalar ktot_min = min(kTotal).value();
                scalar ktot_max = max(kTotal).value();
                Info << "Min kTotal is " << ktot_min << "; max is " << ktot_max << endl;
            }
            else
            {
                Info << "\nTotal turbulent kinetic energy field kTotal already exists!" << endl;
            }
        }
        else
        {
            Info << "\nFluctuating velocity field Uprime not found! Not calculating resolved and total turbulent kinetic energy kResolved and kTotal!" << endl;
        }

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
