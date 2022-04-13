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

            Map<bool> faceMap_;
            EdgeMap<vector> edgeMap_;
            EdgeMap<bool> edgeMapb_;
            pointScalarField alpha1P_;
            Map<Pair<vector>> intFaceMap_;

           //Private Member Functions
           void calculateAlphaP();
           void constructEdgeIndexing();
           void DFS(label i);

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

        //Calculate Interface Area Vector
        tmp<volVectorField::Internal> Sp() const;

        //Calculate Interface Offset
        tmp<volScalarField::Internal> lambda() const;

        //Update Method
        void correct();

        //- Disallow default bitwise assignment
        void operator=(const interfaceReconstruct&) = delete;

        };

} // End namespace Foam
   
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
