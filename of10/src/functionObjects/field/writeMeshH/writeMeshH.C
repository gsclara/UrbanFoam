/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "writeMeshH.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeMeshH, 0);
    addToRunTimeSelectionTable(functionObject, writeMeshH, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeMeshH::writeMeshH
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeMeshH::~writeMeshH()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeMeshH::execute()
{
    return true;
}


bool Foam::functionObjects::writeMeshH::write()
{
    volScalarField V
    (
        IOobject
        (
            mesh_.V().name(),
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedScalar(mesh_.V().name(), mesh_.V().dimensions(), 0),
        calculatedFvPatchField<scalar>::typeName
    );
    
    V.ref() = mesh_.V();

    const dimensionedScalar totalVol
    (
    	sum(V.ref())
    );

    const dimensionedScalar nCell
    (
    	V.size()
    );

    const dimensionedScalar h_cell
    (
    	pow((1/nCell)*totalVol.value(),1.0/3.0)
    );

    Log << "    Writing cell-volumes field " << V.name()
        << " to " << time_.timeName() << endl;
    Log << "    Sum from all cell volumes " << totalVol.value() << endl;
    Log << "    Total number of cells " << V.size() << endl;
    Log << "    h = " << h_cell.value() << endl;

    V.write();


    return true;
}


// ************************************************************************* //
