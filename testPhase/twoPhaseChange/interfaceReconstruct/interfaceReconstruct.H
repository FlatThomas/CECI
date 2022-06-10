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
    Foam::interfaceReconstruct

Description
    Class for Interface Reconstruction based off of contour method,
    using depth first search and hash tables for quick edge indexing/addressing
SourceFiles
    interfaceReconstruct.C

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "compressibleTwoPhaseMixture.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "edgeFaceCirculator.H"

#ifndef interfaceReconstruct_H
#define interfaceReconstruct_H

    namespace Foam
    {
        class interfaceReconstruct
        {
            const compressibleTwoPhaseMixture& mixture_;
            const fvMesh &mesh_;
            pointMesh pMesh_;
            volPointInterpolation vpi_;
            pointScalarField alpha1P_;
            volVectorField::Internal Sp_;
            volScalarField::Internal lambda_;
            volVectorField::Internal intNormal_;

            EdgeMap<label> edgeIndex_;
            Map<bool> faceMap_;
            EdgeMap<bool> boundaryEdges_;
            EdgeMap<vector> edgeMap_;
            EdgeMap<bool> edgeMapb_;
            labelListList edgeFaces_;

           //Private Member Functions
           void calculateAlphaP();
           void createEdgeIndexing();
           void createFaceList();
           void DFS(label i);
           void findBoundaryEdges();
           void generateEdgeFaceList();
           void calculateSp();
           void calculatelambda();
           void calculateintNormal();


        public:
       TypeName("interfaceReconstruct");

        //Construct from Mixture
        interfaceReconstruct(
        const compressibleTwoPhaseMixture &mixture);

        interfaceReconstruct(const interfaceProperties&) = delete;

        //- Destructor
        virtual ~interfaceReconstruct()
        {} 

        const volScalarField::Internal &alpha1() const
            {
                return mixture_.alpha1();
            }

        const volScalarField::Internal &alpha2() const
            {
                return mixture_.alpha2();
            }

        const volVectorField::Internal &Sp() const 
            {
                return Sp_;
            }
        const volScalarField::Internal &lambda() const 
            {
                return lambda_;
            }
        const volVectorField::Internal &intNormal() const 
            {
                return intNormal_;
            }

        
        //Pass Data from an interface field to n number of cells
        template<typename Type, class Mesh>
        void pass(DimensionedField<Type, Mesh> &intField, label n) const;

        //Calculate Interface Area Vector
        void resetFaceMap();
        
        //Update Method
        void correct();

        //- Disallow default bitwise assignment
        void operator=(const interfaceReconstruct&) = delete;

        };

} // End namespace Foam
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
