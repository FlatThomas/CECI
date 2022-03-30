/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "forThermo.H"
#include "makeThermo.H"

#include "specie.H"

#include "thermo.H"

// EoS
#include "icoPolynomial.H"

// Thermo
#include "hPolynomialThermo.H"
#include "sensibleEnthalpy.H"

// Transport
#include "polynomialTransport.H"

// psi/rho
#include "rhoThermo.H"
#include "heRhoThermo.H"

// Mixture
#include "pureMixture.H"


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = c92121a83e9f397179ef59322f976a6a75104c3e
    //
    // Unique function name that can be checked if the correct library version
    // has been loaded
    void heRhoThermo_pureMixture_polynomial_hPolynomial_icoPolynomial_specie___sensibleEnthalpy____c92121a83e9f397179ef59322f976a6a75104c3e(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forThermo
    (
        polynomialTransport,
        sensibleEnthalpy,
        hPolynomialThermo,
        icoPolynomial,
        specie,
        makeThermo,
        rhoThermo,
        heRhoThermo,
        pureMixture
    );
}

// ************************************************************************* //

