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
    Foam::twoPhaseChangeModels::tempGrad

Description
   Temp Gradient Based Model Based off Kunkleman's Thesis 

SourceFiles
    tempGrad.C

\*---------------------------------------------------------------------------*/

#ifndef tempGrad_H
#define tempGrad_H

#include "twoPhaseChangeModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace twoPhaseChangeModels
    {

        /*---------------------------------------------------------------------------*\
                                 Class tempGrad
\*---------------------------------------------------------------------------*/

        class tempGrad
            : public twoPhaseChangeModel
        {
            dimensionedScalar T0;
            dimensionedScalar T1;
            dimensionedScalar T2;
            dimensionedScalar T3;

        public:
            //- Runtime type information
            TypeName("tempGrad");

            // Constructors

            //- Construct for mixture
            tempGrad(
                const compressibleTwoPhaseMixture &mixture);

            //- Destructor
            virtual ~tempGrad()
            {
            }

            // Member Functions
            const volScalarField::Internal &alpha1() const
            {
                return mixture_.alpha1();
            }

            const volScalarField::Internal &alpha2() const
            {
                return mixture_.alpha2();
            }

            const volScalarField::Internal &T() const
            {
                const volScalarField::Internal &T =
                    mixture_.alpha1().db().lookupObject<volScalarField>("T");
                return T;
            }

            const volScalarField::Internal &p() const
            {
                const volScalarField::Internal &p =
                    mixture_.alpha1().db().lookupObject<volScalarField>("p");
                return p;
            }

            //Calculate Liquid Gradient 
            Foam::tmp<Foam::volScalarField::Internal> &liquidGradient() const
            {
                const volScalarField::Internal &alpha1_=alpha1();
                const Foam::fvMesh &mesh=alpha1_.mesh();
                const volScalarField::Internal &Temp=T();


                //Find Liquid Cells that are near interface
                //alpha==1 && neighbour(need to figure out how this works..
                const labelUList& own=mesh.owner();
                const labelUList& nei=mesh.neighbour();
                const volScalarField& interface(mixture_.nearInterface().ref());
                tmp<volScalarField::Internal> tNearIntT(
                        volScalarField::Internal::New(
                            "nearIntT",
                            mesh,
                            dimensionedScalar(dimensionSet(0,0,0,1,0,0,0),0)));
                           
                volScalarField::Internal &nearIntT=tNearIntT.ref();
                
                forAll(own, faceI)
                {
                    //Owner is at Interface, check neighbour alpha
                    if(interface[own[faceI]]!=0)
                    {
                        if(alpha1_[nei[faceI]]==1)
                        {
                            nearIntT[nei[faceI]]=Temp[nei[faceI]];                            
                        }

                    }
                    //Neighbour is at Interface, check owner alpha
                    if(interface[nei[faceI]]!=0)
                    {
                        if(alpha1_[own[faceI]]==1)
                        {
                            nearIntT[own[faceI]]=Temp[own[faceI]];
                        }

                    }

                }
 
               //Shortest distance from cell center to interface


               //Calculate Temp Gradient
    
                forAll(nearIntT,cellI)
                {
                    Info<<nearIntT[cellI]<<endl;
                    Info<<alpha1_[cellI]<<endl;
                }
                return tNearIntT;
            }
                
            
            tmp<volScalarField::Internal>
            Tsat(const volScalarField::Internal &P) const
            {

                tmp<volScalarField::Internal> tmpTsat(
                    new volScalarField::Internal(
                        IOobject(
                            "Tsat",
                            P.time().timeName(),
                            P.mesh(),
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE),

                        P.mesh(),
                        dimensionedScalar("Tsat", dimensionSet(0,0,0, 1,0,0,0), 0)));

                volScalarField::Internal &Tsat = tmpTsat.ref();

                    Tsat=T0+T1*P+T2*pow(P,2)+T3*pow(P,3);

                return tmpTsat;
            }

            //- Return the mass condensation and vaporisation rates as a
            //  coefficient to multiply (1 - alphal) for the condensation rate
            //  and a coefficient to multiply  alphal for the vaporisation rate
            virtual Pair<tmp<volScalarField::Internal>> mDotAlphal();

            //- Return the mass condensation and vaporisation rates as coefficients
            //  to multiply (p - pSat)
            virtual Pair<tmp<volScalarField>> mDotP() const;

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

            //- Correct the tempGrad phaseChange model
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
