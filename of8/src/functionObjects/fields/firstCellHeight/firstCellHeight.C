/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "firstCellHeight.H"
#include "momentumTransportModel.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(firstCellHeight, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        firstCellHeight,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::firstCellHeight::writeFileHeader(const label i)
{
    writeHeader(file(), "firstCellHeight ()");

    writeCommented(file(), "Time");
    writeTabbed(file(), "patch");
    writeTabbed(file(), "min");
    writeTabbed(file(), "max");
    writeTabbed(file(), "average");
    file() << endl;
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::firstCellHeight::calcZfirst
(
    const momentumTransportModel& turbModel
)
{
    tmp<volScalarField> tfirstCellHeight
    (
        volScalarField::New
        (
            type(),
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    volScalarField::Boundary& firstCellHeightBf = tfirstCellHeight.ref().boundaryFieldRef();

    volScalarField::Boundary d = nearWallDist(mesh_).y();

    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
    	firstCellHeightBf[patchi] = d[patchi];
    }

    return tfirstCellHeight;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::firstCellHeight::firstCellHeight
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    writeLocalObjects(obr_, log)
{
    read(dict);
    resetName(typeName);
    resetLocalObjectName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::firstCellHeight::~firstCellHeight()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::firstCellHeight::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    return true;
}


bool Foam::functionObjects::firstCellHeight::execute()
{
    if (mesh_.foundObject<momentumTransportModel>(momentumTransportModel::typeName))
    {
        const momentumTransportModel& model = mesh_.lookupObject<momentumTransportModel>
        (
            momentumTransportModel::typeName
        );

        word name(type());

        return store(name, calcZfirst(model));
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find turbulence model in the "
            << "database" << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::firstCellHeight::write()
{
    Log << type() << " " << name() << " write:" << nl;

    writeLocalObjects::write();

    logFiles::write();

    const volScalarField& firstCellHeight =
        mesh_.lookupObject<volScalarField>(type());

    const volScalarField::Boundary& firstCellHeightBf = firstCellHeight.boundaryField();
    const fvPatchList& patches = mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& patch = patches[patchi];

        if (isA<wallFvPatch>(patch))
        {
            const scalarField& firstCellHeightp = firstCellHeightBf[patchi];

            const scalar minZheight = gMin(firstCellHeightp);
            const scalar maxZheight = gMax(firstCellHeightp);
            const scalar avgZheight = gAverage(firstCellHeightp);

            if (Pstream::master())
            {
                Log << "    patch " << patch.name()
                    << " firstCellHeight : min = " << minZheight << ", max = " << maxZheight
                    << ", average = " << avgZheight << nl;

                writeTime(file());
                file()
                    << tab << patch.name()
                    << tab << minZheight
                    << tab << maxZheight
                    << tab << avgZheight
                    << endl;
            }
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
