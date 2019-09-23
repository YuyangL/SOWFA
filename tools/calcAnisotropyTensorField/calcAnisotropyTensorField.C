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
    calcAnisotropyTensorField

Description
    Calculate and create the anisotropy tensor bij field from exisiting Reynolds stress uuPrime2 field.

Source files:
    readFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

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

        Info << "\nRetrieving anisotropy tensor field bij..." << endl;
        volSymmTensorField bij
        (
            IOobject
            (
                "bij",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            uuPrime2/tr(uuPrime2) - 1/3.*I
        );

        if (!IOobject("bij", runTime.timeName(), mesh).headerOk())
        {
            Info << "\nWriting anisotropy tensor field bij..." << endl;
            bij.write();
        }
        else
        {
            Info << "\nAnisotropy tensor field bij already exists!" << endl;
        }

    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
