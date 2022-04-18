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
        addToRunTimeSelectionTable(twoPhaseChangeModel, tempGrad, interfaceReconstruct);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeModels::tempGrad::tempGrad(
    const compressibleTwoPhaseMixture &mixture, const interfaceReconstruct &interfaceR)
    : twoPhaseChangeModel(typeName, mixture),
    _Reconstruct(interfaceR)
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

Foam::tmp<Foam::volScalarField::Internal>
& Foam::twoPhaseChangeModels::tempGrad::liquidGradient() const
{
    //Standard Tmp Initialization
    Info<<"Initializing tmpGrad field "<<endl;
    Foam::tmp<Foam::volScalarField::Internal> tmpGrad(
        new Foam::volScalarField::Internal(
            IOobject(
                "liqGrad",
                alpha1().mesh().time().timeName(),
                alpha1().mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
                alpha1().mesh(),
                dimensionedScalar("liqGrad", dimensionSet(0,-3,0,1,0,0,0), 0)
            )
        );
        volScalarField::Internal &Grad=tmpGrad.ref();

    Info<<"DONE"<<endl<<endl; 

    //Get References
    const volScalarField::Internal &lambda=_Reconstruct.lambda().ref();
    const fvMesh &mesh=alpha1().mesh();
    const volScalarField::Internal &atInterface=mixture_.nearInterface().ref();
    const labelListList &cellCells=mesh.cellCells();
    const volVectorField::Internal &alphaN=mixture_.n().ref();
    const volScalarField::Internal &satTemp=Tsat(p()).ref();
    

    //Other Variable Declarations
    labelList counter(mesh.C().size());
    labelList interfaceCells;
    labelList dataCells;
    
   Info<<"Find Cells Straddling Interface "<<endl;
   forAll(mesh.C(),cellI)
   {
       if(atInterface[cellI]<1 && lambda[cellI]>1e-4 && alpha1()[cellI]==1)
       {
           Info<<"Cell "<<cellI<<"near Interface "<<endl;
           labelList nearbyIntCells;
           //Cell is Straddling Interface, Find neighbouring interface cells
           forAll(cellCells[cellI],cellJ)
           {
               if(atInterface[cellCells[cellI][cellJ]]==1)
               {
                   Info<<"Neighbour "<<cellJ<<" of "<<cellI<<" at interface"<<endl;
                   nearbyIntCells.append(cellCells[cellI][cellJ]);
                   counter[cellCells[cellI][cellJ]]++;
               }
              
           }

           Info<<"calculating tempGrad at ID'd cells"<<endl;
           forAll(nearbyIntCells,cellJ)
           {
               Grad[cellJ]+=mag(((T()[cellI]-satTemp[nearbyIntCells[cellJ]])/lambda[cellI])
               *alphaN[cellI]);
           }
       }
       //Flag Interface Cells for faster looping down the line
       if(atInterface[cellI]==1)
       {
           interfaceCells.append(cellI);
       }
        
   }

   forAll(interfaceCells,cellI)
   {
        labelList nearbyCells;
        bool breaker=0;
            
        //Loop through neighbouring cells
        forAll(cellCells[interfaceCells[cellI]],cellJ)
        {
            //Break Loop if cell has a liquid neighbour and flip switch 
            if(alpha1()[cellCells[cellI][cellJ]]==1)
            {
                breaker=true;
                break;
            }
            //Also see if cell has neighbour at interface
            else if(atInterface[cellCells[cellI][cellJ]]==1)
            {
                nearbyCells.append(cellCells[cellI][cellJ]);
                counter[cellCells[cellI][cellJ]]++;
            }
        } 
        
        //Only Assign values of breaker wasn't flipped
        if(breaker==false)
        {
            forAll(nearbyCells,cellJ)
            {
                //Only assign grad if breaker was not flipped
                Grad[cellI]+=Grad[nearbyCells[cellJ]];
            }
        }
    }
     

    Info<<"calculating average grad value at cell"<<endl;
    forAll(interfaceCells,cellI)
    {
        if(counter[interfaceCells[cellI]]!=0)
        {
            Grad[interfaceCells[cellI]]/=counter[interfaceCells[cellI]];
        }
    }
    
    return tmpGrad;
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
