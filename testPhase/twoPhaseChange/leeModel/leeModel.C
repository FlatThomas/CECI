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

#include "leeModel.H"
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
        defineTypeNameAndDebug(leeModel, 0);
        addToRunTimeSelectionTable(twoPhaseChangeModel, leeModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::leeModel::leeModel(
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
Foam::twoPhaseChangeModels::leeModel::mDotAlphal()
{
    return Pair<tmp<volScalarField::Internal>>(
        tmp<volScalarField::Internal>(nullptr),
        tmp<volScalarField::Internal>(nullptr));
}



    Foam::Pair<Foam::tmp<Foam::volScalarField>>
    Foam::twoPhaseChangeModels::leeModel::mDotP() const
{
        return Pair<tmp<volScalarField>>(
        tmp<volScalarField>(nullptr),
        tmp<volScalarField>(nullptr));
}
   

Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::leeModel::Salpha(
    volScalarField &alpha) const
{
    //Reference to Mesh
    const volScalarField::Internal &alpha1_ = alpha1();
    const volScalarField::Internal &alpha2_ = alpha2();
    const Foam::fvMesh &mesh = alpha1_.mesh();
        //Initialize Rhodot
    tmp<volScalarField::Internal>trhodotl(
        volScalarField::Internal::New(
            "rhodotl",
            mesh,
            dimensionedScalar(dimMass/ dimTime, 0)));
    
    volScalarField::Internal &rhodotl = trhodotl.ref();

   //Set Interface Field 
   volScalarField::Internal interface(alpha1_); 
   forAll(interface, CellI)
                {
                    if (alpha1_[CellI] < .999 && alpha1_[CellI] > 1e-5)
                    {
                        interface[CellI] = alpha1_[CellI];
                    }
                    else
                    {
                        interface[CellI]=0;
                    }
                }

   const volScalarField::Internal &satTemp=Tsat(p()).ref();
   
   forAll(mesh.C(), CellI)
    {
        //Ensure Cell is at Interface
        if (interface[CellI] > 1e-4)
        {
            //Evaporation
            if (T()[CellI] > satTemp[CellI])
            {
                rhodotl[CellI] = alpha1_[CellI] * .1 * rho1()[CellI]
                                 * ((T()[CellI] - satTemp[CellI]) / satTemp[CellI]);
            }
            //Condensation
            else if (T()[CellI] < satTemp[CellI])
            {
                rhodotl[CellI] = -1*alpha2_[CellI] * .1 * rho2()[CellI]
                    * ((satTemp[CellI] - T()[CellI]) / satTemp[CellI]);
            }
        }
    }


    //Divide by volume to get in terms of density
    rhodotl /= mesh.V();
        tmp<volScalarField::Internal> trhodotv(-1 * rhodotl);

    return Pair<tmp<volScalarField::Internal>>(
        trhodotl,
        trhodotv);
}

Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeModels::leeModel::Sp_rgh(
    const volScalarField &rho,
    const volScalarField &gh,
    volScalarField &p_rgh) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(p_rgh, dimVolume / dimTime));
}

Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeModels::leeModel::SU(
    const volScalarField &rho,
    const surfaceScalarField &rhoPhi,
    volVectorField &U) const
{
    return tmp<fvVectorMatrix>(
        new fvVectorMatrix(U, dimMass * dimVelocity / dimTime));
}

void Foam::twoPhaseChangeModels::leeModel::correct()
{
    twoPhaseChangeModel::correct();
}

bool Foam::twoPhaseChangeModels::leeModel::read()
{
    return twoPhaseChangeModel::read();
}

// ************************************************************************* //
