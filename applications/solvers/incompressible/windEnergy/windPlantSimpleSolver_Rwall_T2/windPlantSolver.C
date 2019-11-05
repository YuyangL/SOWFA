/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    buoyantBoussinesqSimpleFoam

Description
    Steadt-state solver for buoyant, turbulent flow of incompressible fluids

    Uses the Boussinesq approximation:
    \f[
        rho_{k} = 1 - beta(T - T_{ref})
    \f]

    where:
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]

    Valid when:
    \f[
        rho_{k} << 1
    \f]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "IFstream.H"
#include "OFstream.H"
#include "wallDist.H"
#include "interpolateSplineXY.H"
#include "interpolateXY.H"
#include "interpolate2D.H"
#include "horizontalAxisWindTurbinesADM.H"
// From bouoyantBoussinesqSimpleFoam.
// Don't need fvOptions as sources and buoyancy forces are manually added
// #include "fvIOoptionsList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createDivSchemeBlendingField.H"
  //#include "createGradP.H"
    // Source terms are initialized as 0 field
    #include "createSourceTerms.H"
    // #include "readTimeControls.H"
    // #include "CourantNo.H"
    // #include "setInitialDeltaT.H"
  //#include "findVerticalCellLevels.H"
  //#include "findVerticalFaceLevels.H"
  //#include "findWindHeight.H"
  //#include "openCellStatisticsFiles.H"
  //#include "openFaceStatisticsFiles.H"
  //#include "openABLStatisticsFiles.H"
    //#include "createAverageFields.H"
    #include "createVorticityQFields.H"
    #include "computeDivergence.H"
    // From buoyantBoussinesqSimpleFoam
    // Don't need fvOptions as sources and buoyancy forces are manually added
    // #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << nl << "\nStarting time loop\n" << endl;

    // Update boundary conditions before starting in case anything needs
    // updating, for example after using mapFields to interpolate initial
    // field.
    U.correctBoundaryConditions();
    phi = linearInterpolate(U) & mesh.Sf();
    T.correctBoundaryConditions();
  //p_rgh.correctBoundaryConditions();
    turbulence->correct();
    // Shear stress at wall is taken care of by Rwall so nut should have no wall function
    Rwall.correctBoundaryConditions();
    qwall.correctBoundaryConditions();

    // Moved outside of run loop since it only needs to be determined during initializaton
    #include "correctSourceTerms.H"

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;
        // // Time step is always 1 in SimpleFoam
        // Info << "Time Step = " << runTime.timeIndex() << endl;

        // #include "readTimeControls.H"
        // #include "CourantNo.H"
        // #include "setDeltaT.H"
        #include "updateDivSchemeBlendingField.H"


        // --- Pressure-velocity SIMPLE loop
        // while (simple.loop())
        // {
        {
            #include "UEqn.H"
            // Predictor of T based on U* solved from p of previous time
            #include "TEqn.H"
            #include "pEqn.H"
            // TEqn above is predictor, this is corrector based on U solved from p of current time
            #include "TEqn.H"
        }

        // // --- Pressure corrector loop
        // while (pimple.correct())
        // {
        //     #include "pEqn.H"
        //     #include "TEqn.H"
        // }

        // --- Compute the velocity flux divergence
        #include "computeDivergence.H"

//          // --- Update the driving pressure gradient
//          #include "correctGradP.H"

        // // --- Update the source terms
        // // Only update source terms once in the beginning as they'll be constant throughout the whole run
        // if (runTime.value() == 0)
        // {
        //     #include "correctSourceTerms.H"
        // }
        // else
        // {
        //     Info << "Momentum Driving Source Term: (" << SourceU[0].x() << " " << SourceU[0].y() << " " << SourceU[0].z() << ") m/s^2" << endl;
        //     Info << "Temperature Driving Source Term: " << SourceT[0] << " K/s" << endl;
        // }


        // // --- Update the turbulence fields
        // if (pimple.turbCorr())
        // {
        turbulence->correct();
        // }

        // --- Update the turbine array
        turbines.update();

        // --- Update the boundary momentum and
        //     temperature flux conditions
        Rwall.correctBoundaryConditions();
        qwall.correctBoundaryConditions();
        // }


        // #include "computeAverageFields.H"
        #include "computeVorticityQ.H"
//      if (runTime.outputTime())
//      {
//          #include "averageFields.H"
//      }

//      #include "statisticsCell.H"
//      #include "statisticsFace.H"
//      #include "statisticsABL.H"

        // Gonna output Reynolds stress R here now that R has been initialized in createFields.H
        if (runTime.outputTime())
        {
            Rij = turbulence->R();
            Rij.write();
        }

        runTime.write();
//      #include "writeGradP.H"

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End" << endl;

    return 0;
}


// ************************************************************************* //
