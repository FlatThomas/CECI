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

void Foam::interfaceReconstruct::createEdgeIndexing()
{
    Info<<"creating edge indexing "<<endl;
    const edgeList &Edges=mesh_.edges();
    //edgeIndex_.resize(Edges.size());

    forAll(Edges,i)
    {
        edgeIndex_.set(Edges[i],i);
        //Info<<"Edge "<<i<<"added to set"<<endl;
    }
}

void Foam::interfaceReconstruct::createFaceList()
{
    Info<<"Creating Interface Faces List"<<endl;
    //Variable Declarations
    const labelUList& own=mesh_.owner();
    const labelUList& nei=mesh_.neighbour();

    const volScalarField::Internal & interface(mixture_.nearInterface().ref());
    
    forAll(own, faceI)
    {
        //Find owner neighbour pairs both at interface
        if(interface[own[faceI]]==1 && interface[nei[faceI]]==1)
        {
            faceMap_.insert(faceI, 0);
            Info<<"Face "<<faceI<<" added to set"<<mesh_.Cf()[faceI]<<endl;

        }
    }

    forAll(mesh_.boundaryMesh(),patchI)
    {
        const polyPatch& patch=mesh_.boundaryMesh()[patchI];
        const labelUList& pFaceCells=patch.faceCells();

        forAll(patch,faceI)
        {
           if(interface[pFaceCells[faceI]]==1)
           {
               faceMap_.insert(patch.start()+faceI,0);
               Info<<"Face "<<patch.start()+faceI<<" from boundary added to set"<<mesh_.Cf()[patch.start()+faceI]<<endl;
           } 
        }
        
    }
    
    faceMap_.shrink();
    //PreAllocate Space on EdgeMaps
    //edgeMapb_.resize(4*faceMap_.size());
    //edgeMap_.resize(2*faceMap_.size());
}

void Foam::interfaceReconstruct::findBoundaryEdges()
{
    Info<<"finding boundary edges"<<endl;

    boundaryEdges_.resize(mesh_.edges().size());

    forAll(mesh_.boundaryMesh(),patchI)
    {
        //Get Boundary Faces on patch
        const faceList & boundaryFaces=mesh_.boundaryMesh()[patchI];
        forAll(boundaryFaces,FaceI)
        {
            const edgeList &E=boundaryFaces[FaceI].edges();
        
            forAll(E,EdgeI)
            {
                boundaryEdges_.set(E[EdgeI],1);
            }
        }
    }
    boundaryEdges_.shrink();
}
void Foam::interfaceReconstruct::generateEdgeFaceList()
{
    Info<<"Generating EdgeFaceList"<<endl;
    const edgeList &Edges=mesh_.edges();
    const faceList &Faces=mesh_.faces();
    edgeFaces_.setSize(Edges.size());

    forAll(Faces,faceI)
    {
       const edgeList &fE=Faces[faceI].edges();
       forAll(fE,EdgeI)
       {
           label edgeName=edgeIndex_[fE[EdgeI]];
           if(edgeFaces_[edgeName]==labelList::null())
           {
               bool isBoundaryEdge=0;
               
               if(boundaryEdges_.found(fE[EdgeI])==true)
               {
                   isBoundaryEdge=1;
               }

                edgeFaceCirculator circ(mesh_,faceI,1,EdgeI,isBoundaryEdge);
                labelList eF;
                for (edgeFaceCirculator iter=circ.begin(); iter != circ.end(); ++iter) 
                {
                   eF.append(iter.faceLabel());
                }

                edgeFaces_[edgeName]=eF;
           }
       }
    }

}

void Foam::interfaceReconstruct::DFS(label i)
{
    //Info<<"starting DFS loop"<<endl;
    //Mark Current Face as visited
    faceMap_.set(i,1);
    const edgeList &fEdges=mesh_.faces()[i].edges();

    //Loop Over Edges in Face
    for(label edgeI=0; edgeI<fEdges.size(); edgeI++)
    {
        
        //Check if Edge is on List
        
        if(edgeMapb_.found(fEdges[edgeI])==false)
        {
            //Info<<"Edge not found, entering inner loop"<<endl;

            //Add to list
            const edge &E1=fEdges[edgeI];
            edgeMapb_.insert(E1,1);

            if((sign(alpha1P_[E1[0]]-.5)+sign(alpha1P_[E1[1]]-.5))==0)    //Check if Face has .5 somewherebetween verts
            {
                //Info<<"Edge meets criteria, calculated interface point"<<endl;
               //Designations to help clean up code 
               const vector &x1=mesh_.points()[E1[0]];
               const vector &x2=mesh_.points()[E1[1]];
               const scalar &a1=alpha1P_[E1[0]];
               const scalar &a2=alpha1P_[E1[1]];
               
               //Calculate Interpolated Value
               vector linInterp=x1+(x2-x1)*((.5-a1)/(a2-a1));
            
               //Info<<"Linear Interpolation Value found at point "<<linInterp[0]<<" "<<linInterp[1]<<endl;
               //Store in hash table
               edgeMap_.insert(E1,linInterp);

            }

            //Jump to an unvisited Face
            //Info<<"finding edge index"<<endl;
            const label &edgeName=edgeIndex_[E1];
            //Info<<"finding connected faces "<<endl;
            const labelList connectedFaces=edgeFaces_[edgeName];
            
            forAll(connectedFaces,faceI)
            {
                if(faceMap_.found(connectedFaces[faceI])==true)
                    {
                        if(faceMap_[connectedFaces[faceI]]==0)
                        {
                            DFS(connectedFaces[faceI]);
                        }
                    }
            }
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
    faceMap_(mesh_.faces().size()),
    alpha1P_(vpi_.interpolate(mixture_.alpha1()))
        
{
    //calculateAlphaP();
    createFaceList();
    createEdgeIndexing();
    findBoundaryEdges();
    generateEdgeFaceList();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //

tmp<volVectorField> Foam::interfaceReconstruct::Sp() const
{
    //Info<<"constructing Sp vol vector field"<<endl;
    //Standard Tmp Initialization
    tmp<volVectorField>tSp(
        volVectorField::New(
            "SpInt",
            mesh_,
            dimensionedVector(dimensionSet(0,2,0,0,0,0,0),vector::zero)));
    volVectorField &Sp=tSp.ref(); 

    const labelUList &own=mesh_.owner();
    const labelUList &nei=mesh_.neighbour();
    const volVectorField &alphaNormal=mixture_.n().ref();
    const faceList &Faces=mesh_.faces();
 
    //Loop through all Elements in faceMap
    forAllConstIter(Map<bool>,faceMap_, iter)
    {
       label F=iter.key();

       List<vector> V;

       forAll(Faces[F].edges(),edgeI)
       {

           if(edgeMap_.found(Faces[F].edges()[edgeI])==true)
           {
                V.append(edgeMap_[Faces[F].edges()[edgeI]]); 
           }
       }

       vector xproduct=0.5*(V[0]^V[1]);
       /*Info<<"Face "<<F<<" Vertice 1 "<<V[0][0]<<" "<<V[0][1]<<" "<<V[0][2]<<endl;
       Info<<"Face "<<F<<" Vertice 2 "<<V[1][0]<<" "<<V[1][1]<<" "<<V[1][2]<<endl;
       Info<<"Face "<<F<<" X-product "<<xproduct[0]<<" "<<xproduct[1]<<" "<<xproduct[2]<<endl;
      */ 

       if(F<own.size())
       {
            //Flip if Wrong Direction
            if((xproduct & alphaNormal[own[F]])<0)
            {
                xproduct=.5*(V[1]^V[0]);
            }

            Sp[own[F]]+=xproduct;
            Sp[nei[F]]-=xproduct;
                    
        }

       else
       {
           forAll(mesh_.boundaryMesh(),patchI)
           {
               if((mesh_.boundaryMesh()[patchI].start()<= F) &&
               (F < mesh_.boundaryMesh()[patchI].start()+mesh_.boundaryMesh()[patchI].size()))
               {
                    const labelUList &pFaceCells=mesh_.boundaryMesh()[patchI].faceCells();
                    //Info<<"Face "<<F<<" in patch "<<patchI<<endl;

                    if((xproduct & alphaNormal[pFaceCells[F-mesh_.boundaryMesh()[patchI].start()]])<0)
                    {
                        xproduct=.5*(V[1]^V[0]);
                    }

                    Sp[pFaceCells[F-mesh_.boundaryMesh()[patchI].start()]]+=xproduct;
               }
           }
       }
    }

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

    //volScalarField::Internal &lambda=tLambda.ref(); 


    return tLambda;   
}
void Foam::interfaceReconstruct::correct()
{
    Info<<"Starting Correction Loop"<<endl;
    DFS(faceMap_.begin().key());
    
}
