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

#include "interfaceReconstruct.H"

namespace Foam
{
    defineTypeNameAndDebug(interfaceReconstruct, 0);
}
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::interfaceReconstruct::calculateAlphaP()
{
    alpha1P_=vpi_.interpolate(mixture_.alpha1());       
}

void Foam::interfaceReconstruct::constructEdgeIndexing()
{
    //Variable Declarations
    const labelUList& own=mesh_.owner();
    const labelUList& nei=mesh_.neighbour();
    const volScalarField::Internal & interface(mixture_.nearInterface().ref());

     
    forAll(own, faceI)
    {
        //Find owner neighbour pairs both at interface
        if(interface[own[faceI]]==1 && interface[nei[faceI]]==1)
        {
            faceMap_.insert(faceI,0);
        }
    }

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceReconstruct::interfaceReconstruct(
    const compressibleTwoPhaseMixture &mixture):
    mixture_(mixture), 
    mesh_(alpha1().mesh()),
    pMesh_(mesh_),
    vpi_(mesh_),
    alpha1P_(vpi_.interpolate(mixture_.alpha1()))
        
{
    //calculateAlphaP();
    constructEdgeIndexing();
}

void Foam::interfaceReconstruct::DFS(label i)
{
    //Mark Current Face as visited
    faceMap_.set(i,1);
    const edgeList &fEdges=mesh_.faces()[i].edges();

    //Loop Over Edges in Face
    for(label edgeI=0; edgeI<fEdges.size(); edgeI++)
    {
        
        //Check if Edge is on List
        if(edgeMapb_.found(fEdges[edgeI])==false)
        {
            //Add to list
            const edge &E1=fEdges[edgeI];
            edgeMapb_.insert(E1,1);
            
            //Circulator Instantation
            edgeFaceCirculator circ(mesh_,i,true,edgeI,false);
            labelList eFaces;

            //Generate FaceList
            for (edgeFaceCirculator iter=circ.begin(); iter != circ.end(); ++iter) 
            {
                //Only add if in set
                if(faceMap_.found(iter.faceLabel())==true)
                {
                    eFaces.append(iter.faceLabel());
                }
            }
            
            if((sign(alpha1P_[E1[0]]-.5)+sign(alpha1P_[E1[1]]-.5))==0)    //Check if Face has .5 somewherebetween verts
            {
               //Designations to help clean up code 
               const vector &x1=mesh_.points()[E1[0]];
               const vector &x2=mesh_.points()[E1[1]];
               const scalar &a1=alpha1P_[E1[0]];
               const scalar &a2=alpha1P_[E1[1]];
               
               //Calculate Interpolated Value
               vector linInterp=x1+(x2-x1)*((.5-a1)/(a2-a1));
                
                //Store in hash table
                //edgeMap_.insert(E1,linInterp);
                forAll(eFaces,faceJ)
                {
                    Pair<vector> Vs;

                    //If allready Allocated, the old value needs to be stored
                    //in addition to the new value
                    if(intFaceMap_.found(eFaces[faceJ])==true)
                    {
                        Vs.first()=intFaceMap_[eFaces[faceJ]].first();
                        Vs.second()=linInterp;
                        intFaceMap_.set(eFaces[faceJ],Vs);
                    }

                    else
                    {
                        Vs.first()=linInterp;
                        intFaceMap_.set(eFaces[faceJ],Vs);
                    }
                }
            }

            //Jump to an unvisited Face
            forAll(eFaces,faceI)
            {
                if(faceMap_[eFaces[faceI]]==0)    
                    {
                        DFS(faceI);
                    }
            }
        }
    }

   
}

tmp<volVectorField::Internal> Foam::interfaceReconstruct::Sp() const
{
    //Standard Tmp Initialization
    tmp<volVectorField::Internal>tSp(
        volVectorField::Internal::New(
            "SpInt",
            mesh_,
            dimensionedVector(dimensionSet(0,2,0,0,0,0,0),vector::zero)));
    volVectorField::Internal &Sp=tSp.ref(); 

  
  
    //Loop through all Elements in faceMap
    for(Map<Pair<vector>>::const_iterator iter=intFaceMap_.begin(); iter!=intFaceMap_.end(); ++iter)
    {

       const label& F=iter.key();
       const vector& v1=intFaceMap_[F].first();
       const vector& v2=intFaceMap_[F].second();
       
       const labelUList &own=mesh_.owner();
       const labelUList &nei=mesh_.neighbour();
       const volVectorField &alphaNormal=mixture_.n().ref();

       vector xproduct=v1^v2;

       //Flip if Wrong Direction
       if((xproduct & alphaNormal[own[F]])<0)
       {
           xproduct*=-1;
       }


       Sp[own[F]]+=xproduct;
       Sp[nei[F]]-=xproduct;

    }


    Sp*=.5;
    return tSp;   
}

tmp<volScalarField::Internal> Foam::interfaceReconstruct::lambda() const
{
    //Standard Tmp initialization
    tmp<volScalarField::Internal>tLambda(
            volScalarField::Internal::New(
                "lambdaInt",
                mesh_,
                dimensionedScalar(dimensionSet(0,2,0,0,0,0,0),0)));

    volScalarField::Internal &lambda=tLambda.ref(); 


    return tLambda;   
}
void Foam::interfaceReconstruct::correct()
{
    DFS(faceMap_.begin().key());
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
