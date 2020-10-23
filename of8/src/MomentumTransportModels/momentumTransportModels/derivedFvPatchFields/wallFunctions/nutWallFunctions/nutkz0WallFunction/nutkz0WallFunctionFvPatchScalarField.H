/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::nutkz0WallFunctionFvPatchScalarField

Description
    Turbulent viscosity wall-function boundary condition
    for high Reynolds number flows based on near-wall turbulence kinetic energy.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutkz0WallFunction;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutkz0WallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutkz0WallFunctionFvPatchScalarField_H
#define nutkz0WallFunctionFvPatchScalarField_H

#include "nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class nutkz0WallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutkz0WallFunctionFvPatchScalarField
:
    public nutkWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

	//- Surface roughness length
        scalarField z0_;

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("nutkz0WallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutkz0WallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutkz0WallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutkz0WallFunctionFvPatchScalarField
        //  onto a new patch
        nutkz0WallFunctionFvPatchScalarField
        (
            const nutkz0WallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        nutkz0WallFunctionFvPatchScalarField
        (
            const nutkz0WallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutkz0WallFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        nutkz0WallFunctionFvPatchScalarField
        (
            const nutkz0WallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutkz0WallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions
        
        // Return z0
        scalarField& z0()
        {
            return z0_;
        }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );
            //
    // I-O

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
