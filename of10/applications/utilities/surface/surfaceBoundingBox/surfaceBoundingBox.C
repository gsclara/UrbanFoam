/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    surfaceBoundingBox

Description
    Prints the bounding box from a surface, so we can set-up blockMesh and 
    snappyHexMesh accordingly

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "argList.H"
#include "OFstream.H"
#include "OBJstream.H"
#include "SortableList.H"
#include "PatchTools.H"
#include "vtkSurfaceWriter.H"

using namespace Foam;

// Does face use valid vertices?
bool validTri
(
    const bool verbose,
    const triSurface& surf,
    const label facei
)
{
    // Simple check on indices ok.

    const labelledTri& f = surf[facei];

    forAll(f, fp)
    {
        if (f[fp] < 0 || f[fp] >= surf.points().size())
        {
            WarningInFunction
                << "triangle " << facei << " vertices " << f
                << " uses point indices outside point range 0.."
                << surf.points().size()-1 << endl;
            return false;
        }
    }

    if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
    {
        WarningInFunction
            << "triangle " << facei
            << " uses non-unique vertices " << f
            << " coords:" << f.points(surf.points())
            << endl;
        return false;
    }

    // duplicate triangle check

    const labelList& fFaces = surf.faceFaces()[facei];

    // Check if faceNeighbours use same points as this face.
    // Note: discards normal information - sides of baffle are merged.
    forAll(fFaces, i)
    {
        label nbrFacei = fFaces[i];

        if (nbrFacei <= facei)
        {
            // lower numbered faces already checked
            continue;
        }

        const labelledTri& nbrF = surf[nbrFacei];

        if
        (
            ((f[0] == nbrF[0]) || (f[0] == nbrF[1]) || (f[0] == nbrF[2]))
         && ((f[1] == nbrF[0]) || (f[1] == nbrF[1]) || (f[1] == nbrF[2]))
         && ((f[2] == nbrF[0]) || (f[2] == nbrF[1]) || (f[2] == nbrF[2]))
        )
        {
            WarningInFunction
                << "triangle " << facei << " vertices " << f
                << " has the same vertices as triangle " << nbrFacei
                << " vertices " << nbrF
                << " coords:" << f.points(surf.points())
                << endl;

            return false;
        }
    }
    return true;
}


labelList countBins
(
    const scalar min,
    const scalar max,
    const label nBins,
    const scalarField& vals
)
{
    scalar dist = nBins/(max - min);

    labelList binCount(nBins, 0);

    forAll(vals, i)
    {
        scalar val = vals[i];

        label index = -1;

        if (Foam::mag(val - min) < small)
        {
            index = 0;
        }
        else if (val >= max - small)
        {
            index = nBins - 1;
        }
        else
        {
            index = label((val - min)*dist);

            if ((index < 0) || (index >= nBins))
            {
                WarningInFunction
                    << "value " << val << " at index " << i
                    << " outside range " << min << " .. " << max << endl;

                if (index < 0)
                {
                    index = 0;
                }
                else
                {
                    index = nBins - 1;
                }
            }
        }
        binCount[index]++;
    }

    return binCount;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"

    argList::validArgs.append("surface file");
    argList::addBoolOption
    (
        "checkSelfIntersection",
        "also check for self-intersection"
    );
    argList::addBoolOption
    (
        "splitNonManifold",
        "split surface along non-manifold edges"
        " (default split is fully disconnected)"
    );
    argList::addBoolOption
    (
        "verbose",
        "verbose operation"
    );
    argList::addBoolOption
    (
        "blockMesh",
        "write vertices/blocks for blockMeshDict"
    );

    argList args(argc, argv);

    const fileName surfFileName = args[1];

    Info<< "Reading surface from " << surfFileName << " ..." << nl << endl;


    // Read
    // ~~~~

    triSurface surf(surfFileName);


    Info<< "Statistics:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    // write bounding box corners
    if (args.optionFound("blockMesh"))
    {
        pointField cornerPts(boundBox(surf.points(), false).points());

        Info<< "// blockMeshDict info" << nl << nl;

        Info<< "vertices\n(" << nl;
        forAll(cornerPts, ptI)
        {
            Info<< "    " << cornerPts[ptI] << nl;
        }

        // number of divisions needs adjustment later
        Info<< ");\n" << nl
            << "blocks\n"
            << "(\n"
            << "    hex (0 1 2 3 4 5 6 7) (10 10 10) simpleGrading (1 1 1)\n"
            << ");\n" << nl;

        Info<< "edges\n();" << nl
            << "patches\n();" << endl;

        Info<< nl << "// end blockMeshDict info" << nl << endl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
