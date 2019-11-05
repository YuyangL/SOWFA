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
    calc_dDevRab_db

Description
    Calculate and create mean deviatoric Reynolds stress j-gradient tensor field ddev(Rab)/db. Can detect SOWFA LES and RANS case automatically.

Source files:
    createFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"

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

        #include "createFields.H"

        // If mean deviatoric total Ryenolds stress j-gradient tensor field dDevRab_db has not been calculated already, calculate it
        if (!IOobject("dDevRab_db", runTime.timeName(), mesh).headerOk())
        {
            // If uuPrime2 exists, it's from SOWFA LES and calculate mean deviatoric total Ryenolds stress j-gradient tensor field dDevRab_db
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

                Info << "\nComputing mean deviatoric SGS Reynolds stress tensor field, devRSGS..." << endl;
                volSymmTensorField devRSGS
                (
                    IOobject
                    (
                        "devRSGS",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    nuSGSmean*twoSymm(fvc::grad(UAvg))  // Note R in OpenFOAM is defined as ui'uj', not -ui'uj'
                );

                Info << "\nComputing mean total deviatoric Reynolds stress tensor field, devR..." << endl;
                devR = dev(uuPrime2) + devRSGS; // devRResolved is RResolved - 1/3*2kResolved*I
            }
            // Otherwise it's SOWFA RANS
            else
            {
                // Read turbulent kinetic energy field
                Info<< "Reading Reynolds stress field Rij from saved fields..." << endl;
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

                Info << "\nComputing mean deviatoric Reynolds stress tensor field, devR..." << endl;
                devR = dev(Rij);
            }

            scalar devR_min = min(cmptMin(devR));
            scalar devR_max = max(cmptMax(devR));
            Info << "Min devR is " << devR_min << "; max is " << devR_max << endl;

            // Component assignment
            Info << "\nCalculating each component of ddev(Rab)/db..." << endl;
            // grad(dev(R11)), grad(dev(R12)), grad(dev(R13)), grad(dev(R22)), grad(dev(R23)), grad(dev(R33))
            // Accessing devR.xx() doesn't work for a field, only a cell, thus using ".component(symmTensor::XX)" for field
            volVectorField gradR11 = fvc::grad(devR.component(symmTensor::XX));
            volVectorField gradR12 = fvc::grad(devR.component(symmTensor::XY));
            volVectorField gradR13 = fvc::grad(devR.component(symmTensor::XZ));
            volVectorField gradR22 = fvc::grad(devR.component(symmTensor::YY));
            volVectorField gradR23 = fvc::grad(devR.component(symmTensor::YZ));
            volVectorField gradR33 = fvc::grad(devR.component(symmTensor::ZZ));

            // TODO: direct assignment doesn't work for some reason
            // dDevRab_db.component(tensor::XX) = gradR11.component(vector::X);
            // dDevRab_db.component(tensor::XY) = gradR12.component(vector::Y);
            // dDevRab_db.component(tensor::XZ) = gradR13.component(vector::Z);

            // Thus have to assign cell by cell
            forAll(dDevRab_db,cellI)
            {
                // ddev(R11)/dx, ddev(R12)/dy, ddev(R13)/dz
                dDevRab_db[cellI].xx() = gradR11[cellI].x();
                dDevRab_db[cellI].xy() = gradR12[cellI].y();
                dDevRab_db[cellI].xz() = gradR13[cellI].z();
                // ddev(R21)/dx, ddev(R22)/dy, ddev(R23)/dz
                dDevRab_db[cellI].yx() = gradR12[cellI].x();
                dDevRab_db[cellI].yy() = gradR22[cellI].y();
                dDevRab_db[cellI].yz() = gradR23[cellI].z();
                // ddev(R31)/dx, ddev(R32)/dy, ddev(R33)/dz
                dDevRab_db[cellI].zx() = gradR13[cellI].x();
                dDevRab_db[cellI].zy() = gradR23[cellI].y();
                dDevRab_db[cellI].zz() = gradR33[cellI].z();
            }

            Info << "\nWriting mean deviatoric Reynolds stress j-gradient tensor field, dDevRab_db..." << endl;
            dDevRab_db.write();

            scalar magmin = min(cmptMin(dDevRab_db));
            scalar magmax = max(cmptMax(dDevRab_db));
            Info << "Min dDevRab_db is " << magmin << "; max is " << magmax << endl;
        }
        else
        {
            Info << "\nMean deviatoric Reynolds stress j-gradient field dDevRab_db already exists!" << endl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
