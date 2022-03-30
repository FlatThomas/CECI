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

Foam::twoPhaseChangeModels::leeModel::leeModel
(
    const compressibleTwoPhaseMixture& mixture
)
:
    twoPhaseChangeModel(typeName, mixture),
    tSatCoeff_(lookup("pSatCoeff"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeModels::leeModel::mDotAlphal() const
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

        
            
    return Pair<tmp<volScalarField>>
    (
        interfacetmp,
        volScalarField::null()
    );
}



Foam::tmp<Foam::volScalarField::Internal>
Foam::twoPhaseChangeModels::leeModel::mDotP() 
{
     tmp<volScalarField::Internal> test=Foam::twoPhaseChangeModels::leeModel::interfaceCells(); 
     return test;  
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeModels::leeModel::Salpha
(
    volScalarField& alpha
) const
{
    
    return Pair<tmp<volScalarField::Internal>>
    (
        tmp<volScalarField::Internal>(nullptr),
        tmp<volScalarField::Internal>(nullptr)
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeModels::leeModel::Sp_rgh
(
    const volScalarField& rho,
    const volScalarField& gh,
    volScalarField& p_rgh
) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(p_rgh, dimVolume/dimTime));
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeModels::leeModel::SU
(
    const volScalarField& rho,
    const surfaceScalarField& rhoPhi,
    volVectorField& U
) const
{
    return tmp<fvVectorMatrix>
    (
        new fvVectorMatrix(U, dimMass*dimVelocity/dimTime)
    );
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
