/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::twoPhaseChangeModels::leeModel

Description
    Dummy phase change model for 'no phase change'

SourceFiles
    leeModel.C

\*---------------------------------------------------------------------------*/

#ifndef leeModel_H
#define leeModel_H

#include "twoPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace twoPhaseChangeModels
    {

        /*---------------------------------------------------------------------------*\
                                 Class leeModel
\*---------------------------------------------------------------------------*/

        class leeModel
            : public twoPhaseChangeModel
        {

            List<scalar> tSatCoeff_;

        public:
            //- Runtime type information
            TypeName("leeModel");

            // Constructors

            //- Construct for mixture
            leeModel(
                const compressibleTwoPhaseMixture &mixture);

            //- Destructor
            virtual ~leeModel()
            {
            }

            // Member Functions

            const volScalarField::Internal &alpha1() const
            {
                return mixture_.alpha1();
            }

            dimensionedScalar Tsat(scalar P)
            {
                dimensionedScalar Tsat("Tsat", dimTemperature, 0);
                Foam::label arraySize=tSatCoeff_.size();

                for (int i = 0; i < arraySize; i++)
                {
                    Tsat += tSatCoeff_[i] * pow(P, i);
                }

                return Tsat;
            }

            tmp<volScalarField::Internal> interfaceCells()
            {
                const volScalarField::Internal &voidFraction = alpha1();

                tmp<volScalarField::Internal> interfacetmp
                (
                    new volScalarField::Internal(
                    IOobject
                        (
                        "interface",
                        voidFraction.time().timeName(),
                        voidFraction.mesh()
                        ),
                    voidFraction.mesh()
                    )
                );

                volScalarField::Internal& interface = interfacetmp.ref();

                forAll(voidFraction.mesh().C(), CellI)
                {
                    if (voidFraction[CellI] < 1 && voidFraction[CellI] > 0)
                    {
                        interface[CellI] = voidFraction[CellI];
                    }
                }

                return interfacetmp;
            }

            

            //- Return the mass condensation and vaporisation rates as a
            //  coefficient to multiply (1 - alphal) for the condensation rate
            //  and a coefficient to multiply  alphal for the vaporisation rate
            virtual Pair<tmp<volScalarField>> mDotAlphal() const;

            //- Return the mass condensation and vaporisation rates as coefficients
            //  to multiply (p - pSat)
            virtual tmp<volScalarField::Internal> mDotP();

            //- Return unallocated explicit and implicit sources
            //  for the phase-fraction equation
            virtual Pair<tmp<volScalarField::Internal>> Salpha(
                volScalarField &) const;

            //- Return the zero source matrix
            //  for the p_rgh pressure equation
            virtual tmp<fvScalarMatrix> Sp_rgh(
                const volScalarField &rho,
                const volScalarField &gh,
                volScalarField &p_rgh) const;

            //- Return the zero source matrix for the momentum equation
            virtual tmp<fvVectorMatrix> SU(
                const volScalarField &rho,
                const surfaceScalarField &rhoPhi,
                volVectorField &U) const;

            //- Correct the leeModel phaseChange model
            virtual void correct();

            //- Read the phaseProperties dictionary and update
            virtual bool read();
        };

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace twoPhaseChangeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
