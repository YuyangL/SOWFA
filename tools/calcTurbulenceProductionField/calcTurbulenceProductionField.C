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
    calcTurbulenceProductionField

Description
    Calculate and create the turbulence production field.

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

        // If UAvg exists, it's from SOWFA LES and calculate mean turbulence production
        if (IOobject("UAvg", runTime.timeName(), mesh).headerOk())
        {
            // If mean turbulence production field, GAvg, has not been calculated already, calculate it
            if (!IOobject("GAvg", runTime.timeName(), mesh).headerOk())
            {
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

                Info << "\nRetrieving mean turbulence production field, GAvg..." << endl;
                volScalarField GAvg
                (
                    IOobject
                    (
                        "GAvg",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    2.*nuSGSmean*magSqr(symm(fvc::grad(UAvg)))
                );

                Info << "\nWriting mean turbulence production field, GAvg..." << endl;
                GAvg.write();
            }
            else
            {
                Info << "\nMean turbulence production field, GAvg, already exists!" << endl;
            }

        }
        // Else if U is found, then it's SOWFA RANS
        else //if (IOobject("U", runTime.timeName(), mesh).headerOk())
        {
            Info<< "Reading Reynolds-averaged velocity field, U\n" << endl;
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

            // Read turbulent viscosity field
            Info<< "Reading Reynolds-averaged turbulent viscosity field, nut..." << endl;
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

            if (!IOobject("GAvg", runTime.timeName(), mesh).headerOk())
            {


                mix start time not known





                Info << "\nRetrieving mean turbulence production field, GAvg..." << endl;
                scalar mix_ratio = 0.;
                scalar t = runTime_.value();
                // scalar mix_ratio = 0.;
                // Mixing ratio of LES-RANS bij, capped on request
                if (t > mix_startTime_.value())
                {
                    mix_ratio = min((t - mix_startTime_.value())*mix_ratio_cap_.value()/mix_duration_.value(),
                    mix_ratio_cap_.value());
                }

                volScalarField GAvg
                (
                    IOobject
                    (
                        "GAvg",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    ((1. - mix_ratio)*nut_*twoSymm(fvc::grad(U))
                    + mix_ratio*(-2.*k_*bij_))
                    && fvc::grad(U))
                );

                Info << "\nWriting mean turbulence production field, GAvg..." << endl;
                GAvg.write();
            }
            else
            {
                Info << "\nMean turbulence production field, GAvg, already exists!" << endl;
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
