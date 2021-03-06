/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
\*---------------------------------------------------------------------------*/

#include "kEpsilonABL.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kEpsilonABL, 0);
addToRunTimeSelectionTable(RASModelABL, kEpsilonABL, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEpsilonABL::kEpsilonABL
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModelABL(modelName, U, phi, transport, turbulenceModelName),

    // Default of lookupOrAddToDict is dimensionless
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.03  // 123456789.0
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.21
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    // sigmak can also be user specified, like in OF 4.x
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    // Initialize MOST Phi_m/h coefficients
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            coeffDict_,
            5.0
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            16.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            16.0
        )
    ),
    // U* will be inferred from Rwall symmTensor field if not provided
    Ustar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ustar",
            coeffDict_,
            123456789.0,  // Placeholder
            dimensionSet(0, 1, -1, 0, 0, 0, 0)  // Not sure if this is any useful since only value is needed
        )
    ),
    // zetaRef will be inferred by U* and qs if not provided
    zetaRef_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "zetaRef",
            coeffDict_,
            123456789.0  // Placeholder
        )
    ),
    // // IRef has to be provided
    // IRef_(readScalar(coeffDict_.lookup("IRef"))),
    // // Calculated source term for TKE transport, same dimension as Dk/Dt
    // Sk_
    // (
    //     dimensioned<scalar>::lookupOrAddToDict
    //     (
    //         "Sk",
    //         coeffDict_,
    //         0.0,  // For neutral stability
    //         // kMin_.dimensions()/dimTime
    //     )
    // ),
    // Calculated coefficient for buoyancy source/sink in epsilon transport
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            coeffDict_,
            0.0  // For neutral stability
        )
    ),

    // Read k, epsilon, and nut if they're present, otherwise, auto create them
    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        autoCreateEpsilon("epsilon", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),
    // Source term volScalarField in TKE transport
    Sk_
    (
        IOobject
        (
            "Sk",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,  // NO_READ and NO_WRITE by default
            false  // Don't register object
        ),
        U.mesh(),
        dimensionedScalar("Sk", epsilon_.dimensions(), 1.0)
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    // // [NEUTRAL STABILITY ONLY]
    // // If Cmu not provided, infer from provided IRef and zRef
    // // IRef = Cmu^(-1/4)*kappa*sqrt(2/3)/ln(zRef/z0)
    // // Use .value() since Cmu is "dimensioned" although all zeros
    // if (Cmu_.value() == 123456789.0)
    // {
    //     Cmu_.value() = pow(1.0/(IRef_/kappa_/0.8164965809277*log(zRef_/z0_)), 4.0);
    // }

    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    // Initialize source term for TKE transport, a single value to be multiplied by Sk later to form a volScalarField
    scalar Sk_val = 0.0;

    // Initialize Obukhov length to 0
    dimensionedScalar L_("L_", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0);
    // If stable/unstable stability, process U*, L, and zetaRef
    // For some reason qs_.z() doesn't work
    if (qs_.component(2).value() != 0.0)
    {
        // If U* not provided, calculate U* by reading Rwall field at surface and do
        // U* = (mean(u'w')^2 + mean(v'w')^2)^(1/4)
        if (Ustar_.value() == 123456789.0)
        {
            // Set up access to the internal Reynolds stress field (all 0 except at surface)
            // Access it from U's database
            const volSymmTensorField& Rwall = U.db().lookupObject<volSymmTensorField>("Rwall");
            const volScalarField Rwall_xz(Rwall.component(symmTensor::XZ));
            const volScalarField Rwall_yz(Rwall.component(symmTensor::YZ));
            // Set up access to mesh
            const fvMesh& mesh = Rwall.mesh();
            // Find the surface patch ID
            const label patchID = mesh.boundaryMesh().findPatchID("lower");
            // Surface area field of the whole mesh
            const surfaceScalarField& magSf = mesh.magSf();
            // Calculate total area of the lower boundary patch
            scalar area = gSum(magSf.boundaryField()[patchID]);
            // Calculate sum(Rwall_i*area_i) of the lower boundary patch, for u'w' and v'w' component
            scalar areaRwall_xz = gSum(magSf.boundaryField()[patchID]*Rwall_xz.boundaryField()[patchID]);
            scalar areaRwall_yz = gSum(magSf.boundaryField()[patchID]*Rwall_yz.boundaryField()[patchID]);
            // Calculate area weighted average of Rwall = sum(Rwall_i*area_i)/sum(area)
            scalar Rwall_xz_avg = areaRwall_xz/area;
            scalar Rwall_yz_avg = areaRwall_yz/area;
            // Calculate U*
            Ustar_.value() = pow(sqr(Rwall_xz_avg) + sqr(Rwall_yz_avg), 0.25);
            Info << " U* is not provided and is inferred from Rwall symmTensor field as " << Ustar_.value() << " m/s" << endl;
        }

        // Obukhov length
        // volScalarField L_ = pow(Ustar_, 3.0)/kappa_/G_buoyant;
        // Note that g is negative itself
        L_.value() = pow(Ustar_.value(), 3.0)*TRef_.value()/kappa_.value()/(g_.value()&qs_.value());
        // If zetaRef not provided, infer zetaRef from U*, L, zRef and qs
        if (zetaRef_.value() == 123456789.0)
        {
            zetaRef_ = zRef_/L_;
            Info << " zetaRef is not provided and is inferred from U*, L, zRef, and qs as " << zetaRef_.value() << ", bounded to (-2, 1)" << endl;
            // Bound zetaRef to (-2, 1) as per Van der Laan (2016)
            zetaRef_.value() = max(-2.0, zetaRef_.value());
            zetaRef_.value() = min(1.0, zetaRef_.value());
        }
    }
    // Else if neutral stability
    else
    {
        zetaRef_.value() = 0.0;
    }

    // Monin-Obukhov stability parameter
    // volScalarField zeta_ = -kappa_*z_&&(1.0/pow(Ustar_, 3.0))&&(1.0/G_buoyant);

    // Calculate TKE transport source term Sk_ and epsilon transport buoyancy source/sink coefficient C3
    // If unstable ABL
    if (zetaRef_.value() < 0.0)
    {
        scalar phiM_ = pow(1 - gamma1_.value()*zetaRef_.value(), -0.25);
        scalar phiH_ = Prt_.value()*pow(1 - gamma2_.value()*zetaRef_.value(), -0.5);
        scalar phiEps_ = 1.0 - zetaRef_.value();
        scalar fst_ = (2.0 - zetaRef_.value()) + gamma1_.value()/2.0*(1.0 - 12.0*zetaRef_.value() + 7.0*sqr(zetaRef_.value())) - sqr(gamma1_.value())/16.0*zetaRef_.value()*(3.0 - 54.0*zetaRef_.value() + 35.0*sqr(zetaRef_.value()));
        scalar feps_ = pow(phiM_, 2.5)*(1.0 - 0.75*gamma1_.value()*zetaRef_.value());
        // Turbulence intensity at reference height
        // volScalarField IRef_ = sqrt(2.0/3.0*k_)/URef_;
        // scalar Ustar_ = URef_*IRef_*Cmu_*pow(phiEps_, 0.25)*pow(phiM_, 0.25);
        // Constant for TKE transport through diffusion
        scalar Ckappa_ = sqr(kappa_.value())/sigmak_.value()/pow(Cmu_.value(), 0.5);

        // Value of Sk_, to be multiplied by Sk to form a volScalarField
        Sk_val = pow(Ustar_.value(), 3.0)/kappa_.value()/L_.value()*(1.0/zetaRef_.value()*(phiM_ - phiEps_) - phiH_/Prt_.value()/phiM_ - Ckappa_/4.0*pow(phiM_, 6.5)*pow(phiEps_, -1.5)*fst_);
        // Coefficient to regularize the buoyancy source/sink G_buoyant in epsilon transport
        C3_.value() = Prt_.value()/zetaRef_.value()*phiM_/phiH_*(C1_.value()*phiM_ - C2_.value()*phiEps_ + (C2_.value() - C1_.value())*pow(phiEps_, -0.5)*feps_);
    }
    // Else if stable ABL
    else if (zetaRef_.value() > 0.0)
    {
        scalar phiM_ = 1 + beta_.value()*zetaRef_.value();
        scalar phiH_ = Prt_.value() + beta_.value()*zetaRef_.value();
        scalar phiEps_ = phiM_ - zetaRef_.value();
        scalar fst_ = (2.0 - zetaRef_.value()) - 2.0*beta_.value()*zetaRef_.value()*(1.0 - 2.0*zetaRef_.value() + 2.0*beta_.value()*zetaRef_.value());
        scalar feps_ = pow(phiM_, -2.5)*(2.0*phiM_ - 1.0);
        // scalar Ustar_ = URef_*IRef_*Cmu_.value()*pow(pihEps_, 0.25)*pow(phiM_, 0.25);
        scalar Ckappa_ = sqr(kappa_.value())/sigmak_.value()/pow(Cmu_.value(), 0.5);

        Sk_val = pow(Ustar_.value(), 3.0)/kappa_.value()/L_.value()*(1 - phiH_/Prt_.value()/phiM_ - Ckappa_/4.0*pow(phiM_, -3.5)*pow(phiEps_, -1.5)*fst_);
        C3_.value() = Prt_.value()/zetaRef_.value()*phiM_/phiH_*(C1_.value()*phiM_ - C2_.value()*phiEps_ + (C2_.value() - C1_.value())*pow(phiEps_, -0.5)*feps_);
    }

    // Convert a scalar to a volScalarField for TKE transport
    Sk_ *= Sk_val;

    printCoeffs();
    Info << "Source term Sk for TKE transport is " << Sk_val << " m^2/s^3" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Reynolds stress tensor and write it to file
// R = 2/3*k*I - nut*(u_i,j + u_j,i)
// twoSymm returns twice the symmetric part of a tensor, in this case, the tensor is u_i,j, and twoSymm becomes (u_i,j + u_j,i)
tmp<volSymmTensorField> kEpsilonABL::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


// Effective stress tensor incl. laminar stresses
// devReff = -(nu + nut)*(u_i,j + u_j,i - 2/3*u_i,i)
tmp<volSymmTensorField> kEpsilonABL::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                // Why dev"Rho"Reff?
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


// The source term for the incompressible momentum equation
// divDevReff(u_i) = -div((nu + nut)u_i,j) - div((nu + nut)(u_i,j - 1/3*u_i,i).T)
// = -div(nuEff*(u_i,j + (u_i,j).T)), -1/3*u_i,i = 0 due to continuity
// = -div(nuEff*2Sij)
// Laplacian of u_i results in 3x1 vector, i.e. same rank as u_i
tmp<fvVectorMatrix> kEpsilonABL::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


// The source term for the compressible momentum equation
// divDevRhoReff(rho, u_i) = rho*divDevRhoReff
// muEff = rho*(nu + nut)
tmp<fvVectorMatrix> kEpsilonABL::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
)
const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool kEpsilonABL::read()
{
    if (RASModelABL::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        // Read sigmak from dictionary
        sigmak_.readIfPresent(coeffDict());
        // The followings are MOST k-epsilon model related coefficients. Does nothing for neutral stability
        beta_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        Ustar_.readIfPresent(coeffDict());
        zetaRef_.readIfPresent(coeffDict());
        C3_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kEpsilonABL::correct()
{
    RASModelABL::correct();

    Info << "Source term Sk for TKE transport is " << max(Sk_) << " m^2/s^3" << endl;

    if (!turbulence_)
    {
        return;
    }

    // Production term P in k and epsilon transport equation, ABL buoyancy term included
    // Gk = 2nut*sum(Sij)^2
    // Gb = 1/TRef*g*nut/Prt*grad(T), note that g is negative itself
    // symm(u_i,j) returns Sij
    // magSqr of a matrix is the Forbenius norm^2 i.e. sum(elem^2)?
    // [LIMITATION] Prt is constant so not ideal for stable atmospheric stability
    volScalarField G(GName(), nut_*2*magSqr(symm(fvc::grad(U_))));
    volScalarField G_buoyant = (1.0/TRef_)*g_&(nut_/Prt_*fvc::grad(T_));
    // G += G_buoyant;
    // // Richardson number proposed by Rodi
    // volScalarField Rf = -G&&(1.0/(G + G_buoyant));

    // // Monin-Obukhov stability parameter
    // volScalarField z_ = mesh.C().component(vector::Z);
    // volScalarField zeta_ = -kappa_*z_&&(1.0/pow(Ustar_, 3.0))&&(1.0/G_buoyant);

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation rate transport equation with buoyancy source/sink G_buoyant regularized by C3
    // div(phi, epsilon) = nabla*(phi*epsilon), phi is velocity flux, comes from Gauss theorem
    // Sp(a, b) is simply ab
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        (C1_*G + C3_*G_buoyant)*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // TKE transport equation with buoyancy source/sink G_buoyant and additional source Sk_
    // Why epsilon_/k_*k_?
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      + G_buoyant
      - fvm::Sp(epsilon_/k_, k_)
      - Sk_
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
    nut_ = Cmu_*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    // ABL update the turbulent thermal diffusity of T transport
    volScalarField& kappat_ = const_cast<volScalarField&>(U().db().lookupObject<volScalarField>(kappatName_));
    kappat_ = nut_/Prt_;

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
