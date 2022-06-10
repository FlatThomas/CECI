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
#include "tempGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam {
namespace twoPhaseChangeModels {
defineTypeNameAndDebug(tempGrad, 0);
addToRunTimeSelectionTable(twoPhaseChangeModel, tempGrad, interfaceReconstruct);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::tempGrad::tempGrad(
  const compressibleTwoPhaseMixture& mixture,
  const interfaceReconstruct& interfaceR)
  : twoPhaseChangeModel(typeName, mixture)
  , _Reconstruct(interfaceR)
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

// Calculate Sharp Mass Transfer Values Values
Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::tempGrad::mDotAlphal()
{
  return Pair<tmp<volScalarField::Internal>>(
    tmp<volScalarField::Internal>(nullptr),
    tmp<volScalarField::Internal>(nullptr));
}

Foam::tmp<Foam::volScalarField::Internal>
Foam::twoPhaseChangeModels::tempGrad::latentHeat(
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
Foam::twoPhaseChangeModels::tempGrad::Tsat(
  const Foam::volScalarField::Internal& P) const
{
  tmp<volScalarField::Internal> tmpTsat(new volScalarField::Internal(
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

Foam::tmp<Foam::volScalarField::Internal>
Foam::twoPhaseChangeModels::tempGrad::Gradient(
  const Foam::volScalarField::Internal& dint,
  bool isLiquid) const
{
  // Standard Tmp Initialization
  Foam::tmp<Foam::volScalarField::Internal> tmpGrad(
    new Foam::volScalarField::Internal(
      IOobject("liqGrad",
               alpha1().mesh().time().timeName(),
               alpha1().mesh(),
               IOobject::NO_READ,
               IOobject::AUTO_WRITE),
      alpha1().mesh(),
      dimensionedScalar("liqGrad", dimensionSet(0, -3, 0, 1, 0, 0, 0), 0)));
  volScalarField::Internal& Grad = tmpGrad.ref();

  // Get References
  const fvMesh& mesh = alpha1().mesh();
  const labelListList& cellCells = mesh.cellCells();

  tmp<volScalarField> tatInterface = mixture_.nearInterface();
  const volScalarField& atInterface = tatInterface.ref();

  //tmp<volVectorField::Internal> talphaN = _Reconstruct.intNormal();
  //const volVectorField::Internal& alphaN = talphaN.ref();
  const volVectorField::Internal& alphaN = _Reconstruct.intNormal();

  tmp<volScalarField::Internal> tsatTemp = Tsat(p());
  const volScalarField::Internal& satTemp = tsatTemp.ref();

  // Other Variable Declarations
  labelList counter(mesh.C().size(), 0);
  labelList interfaceCells(0);
  scalar cutOff=1;
  if(isLiquid) cutOff=.999;
  else cutOff=.001;


  forAll(mesh.C(), cellI)
  {

    if (dint[cellI] > 1e-4 && (isLiquid? alpha1()[cellI]>cutOff : alpha1()[cellI]<cutOff)) {
      labelList nearbyIntCells;
      // Cell is Straddling Interface, Find neighbouring interface cells
      forAll(cellCells[cellI], cellJ)
      {
        if (atInterface[cellCells[cellI][cellJ]] == 1) {
          nearbyIntCells.append(cellCells[cellI][cellJ]);
          counter[cellCells[cellI][cellJ]]++;
        }
      }

      // Info << "calculating tempGrad at ID'd cells" << endl;
      forAll(nearbyIntCells, cellJ)
      {
        Grad[nearbyIntCells[cellJ]] +=
          mag(((T()[cellI] - satTemp[nearbyIntCells[cellJ]]) / dint[cellI]) *
              alphaN[cellI]);
        Info << "Cell " << nearbyIntCells[cellJ] << endl;
        Info << "Temp at cell " << T()[cellI] << endl;
       // Info << "Alpha N at cell " << alphaN[cellI] << endl;
        Info <<"Saturation Temperature "<< satTemp[nearbyIntCells[cellJ]] << endl;
        Info << "gradient at cell " << nearbyIntCells[cellJ] << " is " <<
        Grad[nearbyIntCells[cellJ]] << endl; Info << "count at a value of " <<
        counter[nearbyIntCells[cellJ]] << endl; Info << "dint value of " <<
        dint[cellI] << endl; Info << endl<<endl; 
      }
    }
    // Flag Interface Cells for faster looping down the line
    if (atInterface[cellI]==1) {
      interfaceCells.append(cellI);
    }
  }

  // Reduce Memory
  tsatTemp.clear();

  forAll(interfaceCells, cellI)
  {
    labelList nearbyCells;
    bool breaker = 0;
    const label currentCell = interfaceCells[cellI];
    const labelList& cellNeighbours = cellCells[currentCell];

    // Loop through neighbouring cells
    forAll(cellNeighbours, cellJ)
    {
      // Break Loop if cell has a liquid neighbour and flip switch
      if(isLiquid? alpha1()[cellNeighbours[cellJ]]>cutOff : alpha1()[cellNeighbours[cellJ]]<cutOff){
      //if (alpha1()[cellNeighbours[cellJ]]<lower || alpha1()[cellNeighbours[cellJ]]>upper) {
        // Info << "Cell " << cellNeighbours[cellJ] << " flipped breaker! " <<
        // endl;
        breaker = true;
        break;
      }
    }

    if (breaker == false) {
      forAll(cellNeighbours, cellJ)
      {
        if (atInterface[cellNeighbours[cellJ]] == 1) {
          Grad[currentCell] += Grad[cellNeighbours[cellJ]];
          counter[currentCell]++;
        }
      }
    }
  }
  tatInterface.clear();

  forAll(interfaceCells, cellI)
  {
    if (counter[interfaceCells[cellI]] != 0) {
      Grad[interfaceCells[cellI]] /= counter[interfaceCells[cellI]];
      // Info << "counter value " << counter[interfaceCells[cellI]] << endl;
      // Info << "Gradient at cell " << interfaceCells[cellI] << " is " <<
      // Grad[interfaceCells[cellI]] << endl;
    }
  }

  return tmpGrad;
}

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::tempGrad::RhoDotSharp(
  const volScalarField& alpha) const
{

  const fvMesh& mesh = alpha.mesh();
  /*const tmp<volScalarField::Internal> lam = _Reconstruct.lambda();
  tmp<volVectorField::Internal> tintNormal=_Reconstruct.intNormal();
  const volVectorField::Internal& intNormal=tintNormal.ref();
 */ 

  const volScalarField::Internal dint =
    mag((mesh.C() & _Reconstruct.intNormal()) - _Reconstruct.lambda());
  
  
  // Standard Tmp INitialization
  tmp<volScalarField::Internal> trhoS(new volScalarField::Internal(
    IOobject("rhoS",
             mesh.time().timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE),

    mesh,
    dimensionedScalar("rhoS", dimMass / dimTime, 0)));

  volScalarField::Internal& rhoS = trhoS.ref();

  // Rando Refs
  const volVectorField::Internal& intArea = _Reconstruct.Sp();
  forAll(intArea, CellI)
  {
    if(mag(intArea[CellI])>0.0){
      Info<<"Cell "<<CellI<<" intArea "<<intArea[CellI]<<endl;
    }
  }
 
  // Grab Thermal Conductivity Refs
  const volScalarField& kL = mixture_.thermo1().kappa();
  const volScalarField& kV = mixture_.thermo2().kappa();

  // Calc Temp Grad at interface due to liquid and vapor
  tmp<volScalarField::Internal> tlgrad = Gradient(dint, 1);
  volScalarField::Internal lgrad = tlgrad.ref();
  tmp<volScalarField::Internal> tvgrad = Gradient(dint, 0);
  volScalarField::Internal vgrad = tvgrad.ref();

  forAll(mesh.C(), cellI)
  {
    if (lgrad[cellI] > 1e-99|| vgrad[cellI] > 1e-99) {
      rhoS[cellI] = ((-kL[cellI] * lgrad[cellI] + kV[cellI] * vgrad[cellI]) /
                     (2260.0 * 10e3)) *
                    mag(intArea[cellI]);
      Info << "Interface Area " << mag(intArea[cellI]) << endl;
      Info << "Liquid Heat Xfer " << -kL[cellI] * lgrad[cellI] << endl;
      Info << "Vapor Heat Xfer " << kV[cellI] * vgrad[cellI] << endl;
      Info << "Calculated mass transfer " << rhoS[cellI] << endl;
      Info << endl<<endl;
    }
  }

  rhoS /= mesh.V();

  tmp<volScalarField::Internal> trhoS2(-1 * rhoS);

  return Pair<tmp<volScalarField::Internal>>(trhoS, trhoS2);
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::tempGrad::mDotP() const
{

  return Pair<tmp<volScalarField>>(tmp<volScalarField>(nullptr),
                                   tmp<volScalarField>(nullptr));
}

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::tempGrad::Salpha(volScalarField& alpha) const
{

  return Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>(
    tmp<volScalarField::Internal>(
      new volScalarField::Internal(rhoSmear().internalField() / rho() * alpha)),
    tmp<volScalarField::Internal>(new volScalarField::Internal(IOobject(
      "zero",
      alpha1().time().timeName(),
      alpha1().mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE),
      alpha1().mesh(),
      dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0)
    )));
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeModels::tempGrad::Sp_rgh(const volScalarField& rho,
                                             const volScalarField& gh,
                                             volScalarField& p_rgh) const
{
  const fvMesh& mesh = alpha1().mesh();
  /*tmp<volScalarField> tmpvdot(
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
    (rhoSmear() * mesh.V()) * (dimlessOne / rho1() - dimlessOne / rho2());*/
  
  tmp<Foam::fvScalarMatrix> tnewEq(new fvScalarMatrix(p_rgh, dimVolume /dimTime));
  Foam::fvScalarMatrix& newEq=tnewEq.ref();
  const dimensionedScalar dimlessOne =
    dimensionedScalar("dimlessOne", dimless, 1.0);

  newEq.source()=(rhoSmear() * mesh.V()) * (dimlessOne / rho1() - dimlessOne / rho2());
  return tnewEq;

}

Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeModels::tempGrad::SU(const volScalarField& rho,
                                         const surfaceScalarField& rhoPhi,
                                         volVectorField& U) const
{
  return tmp<fvVectorMatrix>(
    new fvVectorMatrix(U, dimMass * dimVelocity / dimTime));
}

void
Foam::twoPhaseChangeModels::tempGrad::correct()
{
  twoPhaseChangeModel::correct();
}

bool
Foam::twoPhaseChangeModels::tempGrad::read()
{
  return twoPhaseChangeModel::read();
}

// ************************************************************************* //
