/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
Class
    Foam::incompressible::RASModels::kEpsilonABL
Group
    grpIcoRASTurbulence
Description
    Standard k-epsilon turbulence model for incompressible flows.
    Reference:
    \verbatim
        Launder, B.E., and Spalding, D.B.,
        "Mathematical Models of Turbulence",
        Academic Press, 1972.
        Launder, B.E., and Spalding, D.B.,
        "Computational Methods for Turbulent Flows",
        Comp. Meth in Appl Mech & Eng'g, Vol 3, 1974, pp 269-289.
    \endverbatim
    The default model coefficients correspond to the following:
    \verbatim
        kEpsilonABLCoeffs
        {
            Cmu         0.03;  // 0.09;
            C1          1.21;  // 1.44;
            C2          1.92;
            sigmaEps    1.3;
            sigmak      1.0;
            beta        5.0;
        }
    \endverbatim
SourceFiles
    kEpsilonABL.C
\*---------------------------------------------------------------------------*/

#ifndef kEpsilonABL_H
#define kEpsilonABL_H

#include "RASModelABL.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kEpsilonABL Declaration
\*---------------------------------------------------------------------------*/

class kEpsilonABL
:
    public RASModelABL
{

protected:

    // Protected data

        // Standard model coefficients
            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaEps_;
            // Added sigmak
            dimensionedScalar sigmak_;

        // MOST Phi_m/h coefficients
        // beta is coefficient of zeta for stable Phi_m/h
        // gamma1 is coefficient of zeta for unstable Phi_m
        // gamma2 is coefficient of zeta for unstable Phi_h
            dimensionedScalar beta_;
            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

        // Provided k-epsilon MOST model coefficients
        // U* is friction speed, if not provided, U* will be inferred from Rwall symmTensor field
        // zetaRef is reference Monin-Obukhov stability parameter at reference height
        // IRef is reference turbulence intensity at reference height
            dimensionedScalar Ustar_;
            dimensionedScalar zetaRef_;
            // scalar Iref_;

        // Calculated k-epsilon MOST model coefficients,
        // C3 the coefficient of buoyancy source/sink in epsilon transport
            dimensionedScalar C3_;


        // Fields
        // Sk the source term in TKE transport
            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;
            volScalarField Sk_;


public:

    //- Runtime type information
    TypeName("kEpsilonABL");

    // Constructors

        //- Construct from components
        kEpsilonABL
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~kEpsilonABL()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        // According Van der Laan (2014), DkEff = (nut/sigmak + nu)
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/sigmak_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the source term in TKE transport
        virtual volScalarField Sk() const
        {
            return Sk_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
