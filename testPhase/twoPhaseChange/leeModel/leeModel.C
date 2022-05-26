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

#include "addToRunTimeSelectionTable.H"
#include "fvScalarMatrix.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "leeModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam {
namespace twoPhaseChangeModels {
defineTypeNameAndDebug(leeModel, 0);
addToRunTimeSelectionTable(twoPhaseChangeModel, leeModel, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::leeModel::leeModel(
  const compressibleTwoPhaseMixture& mixture)
  : twoPhaseChangeModel(typeName, mixture)
{
  twoPhaseChangeModelCoeffs_.lookup("T0") >> T0;
  twoPhaseChangeModelCoeffs_.lookup("T1") >> T1;
  twoPhaseChangeModelCoeffs_.lookup("T2") >> T2;
  twoPhaseChangeModelCoeffs_.lookup("T3") >> T3;

  twoPhaseChangeModelCoeffs_.lookup("H0") >> H0;
  twoPhaseChangeModelCoeffs_.lookup("H1") >> H1;
  twoPhaseChangeModelCoeffs_.lookup("H2") >> H2;
  twoPhaseChangeModelCoeffs_.lookup("H3") >> H3;
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::leeModel::RhoDotSharp(
  const volScalarField& alpha) const
{
  // Reference to Mesh
  const volScalarField::Internal& alpha1_ = alpha1();
  const volScalarField::Internal& alpha2_ = alpha2();
  const Foam::fvMesh& mesh = alpha1_.mesh();
  // Initialize Rhodot
  tmp<volScalarField::Internal> trhodotl(volScalarField::Internal::New(
    "rhodotl", mesh, dimensionedScalar(dimMass / dimTime, 0)));

  volScalarField::Internal& rhodotl = trhodotl.ref();

  // Set Interface Field
  volScalarField::Internal interface(alpha1_);
  forAll(interface, CellI)
  {
    if (alpha1_[CellI] < .999 && alpha1_[CellI] > 1e-5) {
      interface[CellI] = alpha1_[CellI];
    } else {
      interface[CellI] = 0;
    }
  }

  const volScalarField::Internal& satTemp = Tsat(p()).ref();

  forAll(mesh.C(), CellI)
  {
    // Ensure Cell is at Interface
    if (interface[CellI] > 1e-4) {
      // Evaporation
      if (T()[CellI] > satTemp[CellI]) {
        rhodotl[CellI] = alpha1_[CellI] * .1 * rho1()[CellI] *
                         ((T()[CellI] - satTemp[CellI]) / satTemp[CellI]);
      }
      // Condensation
      else if (T()[CellI] < satTemp[CellI]) {
        rhodotl[CellI] = -1 * alpha2_[CellI] * .1 * rho2()[CellI] *
                         ((satTemp[CellI] - T()[CellI]) / satTemp[CellI]);
      }
    }
  }

  // Divide by volume to get in terms of density
  rhodotl /= mesh.V();
  tmp<volScalarField::Internal> trhodotv(-1 * rhodotl);

  return Pair<tmp<volScalarField::Internal>>(trhodotl, trhodotv);
}

Foam::tmp<Foam::volScalarField::Internal>
Foam::twoPhaseChangeModels::leeModel::latentHeat(
  const Foam::volScalarField::Internal& T) const
{

  tmp<volScalarField::Internal> tmplatentHeat(new volScalarField::Internal(
    IOobject("latentHeat",
             alpha1().time().timeName(),
             alpha1().mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE),
    alpha1().mesh(),
    dimensionedScalar("latentHeat", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0)));

  volScalarField::Internal& latentHeat = tmplatentHeat.ref();
  latentHeat = H0 + H1 * T + H2 * pow(T, 2) + H3 * pow(T, 3);
  return tmplatentHeat;
}
Foam::tmp<Foam::volScalarField::Internal>
Foam::twoPhaseChangeModels::leeModel::Tsat(
  const Foam::volScalarField::Internal& P) const
{

  Foam::tmp<Foam::volScalarField::Internal> tmpTsat(
    new Foam::volScalarField::Internal(
      IOobject("Tsat",
               P.time().timeName(),
               P.mesh(),
               IOobject::NO_READ,
               IOobject::AUTO_WRITE),

      P.mesh(),
      dimensionedScalar("Tsat", dimensionSet(0, 0, 0, 1, 0, 0, 0), 0)));

  volScalarField::Internal& Tsat = tmpTsat.ref();

  Tsat = T0 + T1 * P + T2 * pow(P, 2) + T3 * pow(P, 3);

  return tmpTsat;
}
// Calculate Sharp Mass Transfer Values Values
Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::leeModel::mDotAlphal()
{
  return Pair<tmp<volScalarField::Internal>>(
    tmp<volScalarField::Internal>(nullptr),
    tmp<volScalarField::Internal>(nullptr));
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::leeModel::mDotP() const
{
  return Pair<tmp<volScalarField>>(tmp<volScalarField>(nullptr),
                                   tmp<volScalarField>(nullptr));
}

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::leeModel::Salpha(volScalarField& alpha) const
{

  return Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>(
    tmp<volScalarField::Internal>(
      new volScalarField::Internal(rhoSmear().internalField() / rho() * alpha)),
    tmp<volScalarField::Internal>(nullptr));
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeModels::leeModel::Sp_rgh(const volScalarField& rho,
                                             const volScalarField& gh,
                                             volScalarField& p_rgh) const
{
  const fvMesh& mesh = alpha1().mesh();
  tmp<volScalarField> tmpvdot(
    new volScalarField(IOobject("vdot",
                                alpha1().time().timeName(),
                                alpha1().mesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE),
                       alpha1().mesh(),
                       dimensionedScalar("vdot", dimVolume / dimTime, 0)));

  const dimensionedScalar dimlessOne =
    dimensionedScalar("dimlessOne", dimless, 1.0);

  volScalarField& vdot = tmpvdot.ref();
  vdot.ref() =
    (rhoSmear() * mesh.V()) / (dimlessOne / rho1() - dimlessOne / rho2());

  return tmp<fvScalarMatrix>(
    new fvScalarMatrix(tmpvdot.ref(), dimVolume / dimTime));
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeModels::leeModel::SU(const volScalarField& rho,
                                         const surfaceScalarField& rhoPhi,
                                         volVectorField& U) const
{
  return tmp<fvVectorMatrix>(
    new fvVectorMatrix(U, dimMass * dimVelocity / dimTime));
}

void
Foam::twoPhaseChangeModels::leeModel::correct()
{
  twoPhaseChangeModel::correct();
}

bool
Foam::twoPhaseChangeModels::leeModel::read()
{
  return twoPhaseChangeModel::read();
}

// ************************************************************************* //
