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

\*---------------------------------------------------------------------------*/

#include "tempGrad.H"
#include "fvScalarMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace twoPhaseChangeModels
    {
        defineTypeNameAndDebug(tempGrad, 0);
        addToRunTimeSelectionTable(twoPhaseChangeModel, tempGrad, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::tempGrad::tempGrad(
    const compressibleTwoPhaseMixture &mixture)
    : twoPhaseChangeModel(typeName, mixture)
{
    twoPhaseChangeModelCoeffs_.lookup("T0")>>T0;
    twoPhaseChangeModelCoeffs_.lookup("T1")>>T1;
    twoPhaseChangeModelCoeffs_.lookup("T2")>>T2;
    twoPhaseChangeModelCoeffs_.lookup("T3")>>T3;
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//Calculate Sharp Mass Transfer Values Values
Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::tempGrad::mDotAlphal()
{
    return Pair<tmp<volScalarField::Internal>>(
        tmp<volScalarField::Internal>(nullptr),
        tmp<volScalarField::Internal>(nullptr));
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::tempGrad::mDotP() const
{
        return Pair<tmp<volScalarField>>(
        tmp<volScalarField>(nullptr),
        tmp<volScalarField>(nullptr));
}
   


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::tempGrad::Salpha(
    volScalarField &alpha) const
{
      return Pair<tmp<volScalarField::Internal>>
          (
           tmp<volScalarField::Internal>(nullptr),
           tmp<volScalarField::Internal>(nullptr)
          );
        
        
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeModels::tempGrad::Sp_rgh(
    const volScalarField &rho,
    const volScalarField &gh,
          volScalarField &p_rgh) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(p_rgh, dimVolume / dimTime));
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeModels::tempGrad::SU(
    const volScalarField &rho,
    const surfaceScalarField &rhoPhi,
          volVectorField &U) const
{
    return tmp<fvVectorMatrix>(
        new fvVectorMatrix(U, dimMass * dimVelocity / dimTime));
}

void Foam::twoPhaseChangeModels::tempGrad::correct()
{
    twoPhaseChangeModel::correct();
}

bool Foam::twoPhaseChangeModels::tempGrad::read()
{
    return twoPhaseChangeModel::read();
}

// ************************************************************************* //
