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

namespace Foam {
defineTypeNameAndDebug(interfaceReconstruct, 0);
}
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void
Foam::interfaceReconstruct::calculateAlphaP()
{
  alpha1P_ = vpi_.interpolate(mixture_.alpha1());
}

void
Foam::interfaceReconstruct::createEdgeIndexing()
{
  Info << "creating edge indexing " << endl;
  const edgeList& Edges = mesh_.edges();
  // edgeIndex_.resize(Edges.size());

  forAll(Edges, i)
  {
    edgeIndex_.set(Edges[i], i);
    // Info<<"Edge "<<i<<"added to set"<<endl;
  }
}

void
Foam::interfaceReconstruct::createFaceList()
{
  // Variable Declarations
  const labelUList& own = mesh_.owner();
  const labelUList& nei = mesh_.neighbour();

  const volScalarField::Internal& interface=mixture_.nearInterface().ref();
  faceMap_.clear();

  forAll(own, faceI)
  {
    // Find owner neighbour pairs both at interface
    if (interface[own[faceI]] == 1 && interface[nei[faceI]] == 1) {
      faceMap_.set(faceI, 0);
    }
  }

  forAll(mesh_.boundaryMesh(), patchI)
  {
    const polyPatch& patch = mesh_.boundaryMesh()[patchI];
    const labelUList& pFaceCells = patch.faceCells();

    forAll(patch, faceI)
    {
      if (interface[pFaceCells[faceI]] == 1) {
        faceMap_.set(patch.start() + faceI, 0);
         }
    }
  }

  faceMap_.shrink();
  // PreAllocate Space on EdgeMaps
  // edgeMapb_.resize(4*faceMap_.size());
  // edgeMap_.resize(2*faceMap_.size());
}

void
Foam::interfaceReconstruct::findBoundaryEdges()
{

  boundaryEdges_.resize(mesh_.edges().size());

  forAll(mesh_.boundaryMesh(), patchI)
  {
    // Get Boundary Faces on patch
    const faceList& boundaryFaces = mesh_.boundaryMesh()[patchI];
    forAll(boundaryFaces, FaceI)
    {
      const edgeList& E = boundaryFaces[FaceI].edges();

      forAll(E, EdgeI) { boundaryEdges_.set(E[EdgeI], 1); }
    }
  }
  boundaryEdges_.shrink();
}
void
Foam::interfaceReconstruct::generateEdgeFaceList()
{
  const edgeList& Edges = mesh_.edges();
  const faceList& Faces = mesh_.faces();
  edgeFaces_.setSize(Edges.size());

  forAll(Faces, faceI)
  {
    const edgeList& fE = Faces[faceI].edges();
    forAll(fE, EdgeI)
    {
      label edgeName = edgeIndex_[fE[EdgeI]];
      if (edgeFaces_[edgeName] == labelList::null()) {
        bool isBoundaryEdge = 0;

        if (boundaryEdges_.found(fE[EdgeI]) == true) {
          isBoundaryEdge = 1;
        }

        edgeFaceCirculator circ(mesh_, faceI, 1, EdgeI, isBoundaryEdge);
        labelList eF;
        for (edgeFaceCirculator iter = circ.begin(); iter != circ.end();
             ++iter) {
          eF.append(iter.faceLabel());
        }

        edgeFaces_[edgeName] = eF;
      }
    }
  }
}

void
Foam::interfaceReconstruct::DFS(label i)
{
  // Mark Current Face as visited
  faceMap_.set(i, 1);
  const edgeList& fEdges = mesh_.faces()[i].edges();
  // Loop Over Edges in Face
  for (label edgeI = 0; edgeI < fEdges.size(); edgeI++) {

    // Check if Edge is on List

    if (edgeMapb_.found(fEdges[edgeI]) == false) {

      // Add to list
      const edge& E1 = fEdges[edgeI];
      edgeMapb_.insert(E1, 1);

      if ((sign(alpha1P_[E1[0]] - .5) + sign(alpha1P_[E1[1]] - .5)) == 0 ||
            mag(alpha1P_[E1[0]]-.5)<1e-4 || mag(alpha1P_[E1[1]]-.5)<1e-4) // Check if Face has .5 somewherebetween verts
      {
        // Designations to help clean up code
        const vector& x1 = mesh_.points()[E1[0]];
        const vector& x2 = mesh_.points()[E1[1]];
        const scalar& a1 = alpha1P_[E1[0]];
        const scalar& a2 = alpha1P_[E1[1]];
                // Calculate Interpolated Value
        vector linInterp = x1 + (x2 - x1) * ((.5 - a1) / (a2 - a1));

         //Info<<"Linear Interpolation Value found at point "<<linInterp<<endl;
        edgeMap_.insert(E1, linInterp);
      }
      

      // Jump to an unvisited Face
      // Info<<"finding edge index"<<endl;
      const label& edgeName = edgeIndex_[E1];
      // Info<<"finding connected faces "<<endl;
      const labelList connectedFaces = edgeFaces_[edgeName];

      forAll(connectedFaces, faceI)
      {
        if (faceMap_.found(connectedFaces[faceI]) == true) {
          if (faceMap_[connectedFaces[faceI]] == 0) {
            DFS(connectedFaces[faceI]);
          }
        }
      }
    }
  }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceReconstruct::interfaceReconstruct(
  const compressibleTwoPhaseMixture& mixture)
  : mixture_(mixture)
  , mesh_(alpha1().mesh())
  , pMesh_(mesh_)
  , vpi_(mesh_)
  , faceMap_(mesh_.faces().size())
  , alpha1P_(vpi_.interpolate(mixture_.alpha1()))
  , Sp_
  (
      IOobject
      (
          "Sp",
          mixture_.alpha1().time().timeName(),
          mixture.alpha1().mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mixture.alpha1().mesh(),
      dimensionedVector(dimArea, vector::zero)
  ),

  lambda_
  (
      IOobject
      (
          "lambda",
          mixture_.alpha1().time().timeName(),
          mixture.alpha1().mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mixture.alpha1().mesh(),
      dimensionedScalar(dimLength, 0)
  ),

  intNormal_
  (
      IOobject
      (
          "intNormal",
          mixture_.alpha1().time().timeName(),
          mixture.alpha1().mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      mixture.alpha1().mesh(),
      dimensionedVector(dimless, vector::zero)
  )

{
  // calculateAlphaP();
  createFaceList();
  createEdgeIndexing();
  findBoundaryEdges();
  generateEdgeFaceList();
  DFS(faceMap_.begin().key());
  calculateSp();
  calculateintNormal();
  calculatelambda();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
template<typename Type, class Mesh>
void
Foam::interfaceReconstruct::pass(DimensionedField<Type, Mesh>& intField,
                                 label n) const
{
  const labelListList& cellCells = mesh_.cellCells();
  const volScalarField& atInterface = mixture_.nearInterface().ref();

  labelList intCells(0);
  labelList queueCells(0);
  labelList counter(mesh_.C().size(), 0);
  labelList visited(mesh_.C().size(), 0);

  // Prelim Search
  forAll(mesh_.C(), cellI)
  {
    if (mag(intField[cellI]) > 1e-6) // Maybe a better way to do this
    {
      intCells.append(cellI);
      visited[cellI] = 1;
    }
  }
  // Repeat this cycle n times
  for (int i = 0; i < n; i++) {

    // Info<<"beginning loop "<<i<<" "<<endl;
    // Loop over interface cells
    forAll(intCells, cellI)
    {
      labelList adjCells =
        cellCells[intCells[cellI]]; // ID cells adjacent to current cell
      // Start Loop over adjacent Cells
      // Info<<"Starting loop over adjacent cells"<<endl;
      // Info<<"..."<<endl;

      forAll(adjCells, cellJ)
      {
        // Look for cells not already visited
        if (visited[adjCells[cellJ]] == 0) {
          counter[adjCells[cellJ]]++; // Counter for average

          // Add to index if first time being visited
          if (counter[adjCells[cellJ]] == 1) {
            queueCells.append(adjCells[cellJ]);
          }

          intField[adjCells[cellJ]] +=
            intField[intCells[cellI]]; // add interface value to adjacent cell
        }
      }
      // Info<<"DONE"<<endl;
    }

    // Clear and Resize Interface Cell for next major loop
    intCells.clear();
    intCells.resize(queueCells.size());

    // Loop through Queued cells and calculate Average
    forAll(queueCells, cellI)
    {

      intField[queueCells[cellI]] /= counter[queueCells[cellI]];

      // Mark queued cells as visited
      visited[queueCells[cellI]] = 1;

      // Move Queued Cells to Interface
      intCells[cellI] = queueCells[cellI];
    }

    queueCells.clear();
  }

  
}

void Foam::interfaceReconstruct::calculateSp() 
{
  // Standard Tmp Initialization
  
  const labelUList& own = mesh_.owner();
  const labelUList& nei = mesh_.neighbour();
  const volVectorField& alphaNormal = mixture_.n().ref();
  const faceList& Faces = mesh_.faces();
  Sp_=zeroField();

  // Loop through all Elements in faceMap
  forAllConstIter(Map<bool>, faceMap_, iter)
  {
    label F = iter.key();

    List<vector> V;

    forAll(Faces[F].edges(), edgeI)
    {
      if (edgeMap_.found(Faces[F].edges()[edgeI]) == true) {

        V.append(edgeMap_[Faces[F].edges()[edgeI]]);
      }
    }

    if(V.empty())
      {
        continue;
      }
    vector xproduct = 0.5 * (V[0] ^ V[1]);

    if (F < own.size()) {
      // Flip if Wrong Direction
      if ((xproduct & alphaNormal[own[F]]) <= 0) {
        xproduct = .5 * (V[1] ^ V[0]);
      }

      Sp_[own[F]] += xproduct;
      Sp_[nei[F]] -= xproduct;
    }

    else {
      forAll(mesh_.boundaryMesh(), patchI)
      {
        if ((mesh_.boundaryMesh()[patchI].start() <= F) &&
            (F < mesh_.boundaryMesh()[patchI].start() +
                   mesh_.boundaryMesh()[patchI].size())) {
          const labelUList& pFaceCells =
            mesh_.boundaryMesh()[patchI].faceCells();
          // Info<<"Face "<<F<<" in patch "<<patchI<<endl;

          if ((xproduct &
               alphaNormal[pFaceCells[F -
                                      mesh_.boundaryMesh()[patchI].start()]]) <=
              0) {
            xproduct = .5 * (V[1] ^ V[0]);
          }
          Sp_[pFaceCells[F - mesh_.boundaryMesh()[patchI].start()]] += xproduct;
        }
      }
    }
  }

  pass(Sp_, 1);

}


void Foam::interfaceReconstruct::calculatelambda() 
{
  // Standard Tmp initialization
  labelList counter(mesh_.nCells(), 0);

  const labelUList& own = mesh_.owner();
  const labelUList& nei = mesh_.neighbour();
  const volVectorField::Internal& alphGrad = intNormal();
  lambda_=zeroField();

  for (HashTable<vector, edge, Hash<edge>>::const_iterator iter =
         edgeMap_.cbegin();
       iter != edgeMap_.end();
       iter++) {
    edge eI = iter.key();
    label edgeLabel = edgeIndex_[eI];
    const labelList& eFaces = edgeFaces_[edgeLabel];

    forAll(eFaces, faceI)
    {
      label F = eFaces[faceI];

      if (eFaces[faceI] < own.size()) {
        lambda_[own[F]] += alphGrad[own[F]] & edgeMap_[iter.key()];
        lambda_[nei[F]] += alphGrad[nei[F]] & edgeMap_[iter.key()];
        counter[own[F]]++;
        counter[nei[F]]++;
      }

      else {
        forAll(mesh_.boundaryMesh(), patchI)
        {
          if ((mesh_.boundaryMesh()[patchI].start() <= F) &&
              (F < mesh_.boundaryMesh()[patchI].start() +
                     mesh_.boundaryMesh()[patchI].size())) {
            const labelUList& pFaceCells =
              mesh_.boundaryMesh()[patchI].faceCells();
            label boundaryCellIndex =
              pFaceCells[F - mesh_.boundaryMesh()[patchI].start()];

            lambda_[boundaryCellIndex] +=
              alphGrad[boundaryCellIndex] & edgeMap_[iter.key()];
            counter[boundaryCellIndex]++;
            
          }
        }
      }
    }
  }

  // Loop over Cells and divide by counter value

  forAll(mesh_.C(), cellI)
  {
    if (counter[cellI] != 0) {
      lambda_[cellI] /= counter[cellI];
    }
  }
  pass(lambda_, 1);

}

 
void Foam::interfaceReconstruct::calculateintNormal() 
{

  forAll(intNormal_, cellI)
  {
    if(mag(Sp()[cellI])>1e-6)
    {
      intNormal_[cellI]=Sp()[cellI]/mag(Sp()[cellI]);
    }
  }
}

void Foam::interfaceReconstruct::correct()
{
  calculateAlphaP();
  createFaceList();
  edgeMap_.clearStorage();
  edgeMapb_.clearStorage();
  DFS(faceMap_.begin().key());
  calculateSp();
  calculateintNormal();
  calculatelambda();
}
