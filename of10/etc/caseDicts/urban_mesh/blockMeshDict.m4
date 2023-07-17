/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     |                                                 |
|   \\  /    A nd           | Copyright (C) 2016 Ehsan Madadi-Kandjani        |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// General macros to create cylinder mesh
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])
define(hex2D0, hex ($1b $2b $3b $4b $1m $2m $3m $4m))
define(hex2D1, hex ($1m $2m $3m $4m $1t $2t $3t $4t))
define(btQuad0, ($1b $2b $1m $2m))
define(btQuad1, ($1m $2m $2t $1t))
define(topQuad, ($1t $4t $3t $2t))
define(bottomQuad, ($1b $2b $3b $4b))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Maximum height in the domain
define(H, 32)

// Center point in X cylinder
define(x0, 0)
// Center point in Y cylinder
define(y0, 0)
// Define bottom height (from the terrain)
define(z0,-5.62434 )

// Inner square side half //todo: rewrite
define(s, 720)
// Inner square side curvature
//define(sc, 720) //todo: try to parameterise
define(sc, calc(s*1.1)) //todo: try to parameterise

// Number of cells at inner square
define(Ns0, 56)
// Number of cells between inner square and circle
define(Ni0, 28)
// Number of cells in the cylinder height
define(Nz0, 4)
define(Nz1, 6)

// Define stretching 
define(stretchZ, 2)
define(stretchXY, 2)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Calculate CFD domain extension following COST732 guidelines
define(Rcity, calc(15*H + 2*s) )
//define(Rcity, 1500 )
// Cylinder radius
define(r, Rcity)

// Second part of cells in inner/outer circle
define(Ns1, Ns0)
define(Ni1, Ni0)

// Height of cylinder
define(z, calc(10*H))
// Base z
define(Zb, z0)
// Mid z
define(Zm, calc(Zb + 2*H))
// Outlet z
define(Zt, calc(Zb + z))

define(vert, (x$1$2 y$1$2 $3))
define(evert, (ex$1$2 ey$1$2 $3))

// 45 degree points angle
define(a0, -45)
define(a1, -135)
define(a2, 135)
define(a3, 45)

// Half of 45 degree points angle
define(ea0, 0)
define(ea1, -90)
define(ea2, 180)
define(ea3, 90)

define(ca0, calc(cos((pi/180)*a0)))
define(ca1, calc(cos((pi/180)*a1)))
define(ca2, calc(cos((pi/180)*a2)))
define(ca3, calc(cos((pi/180)*a3)))

define(sa0, calc(sin((pi/180)*a0)))
define(sa1, calc(sin((pi/180)*a1)))
define(sa2, calc(sin((pi/180)*a2)))
define(sa3, calc(sin((pi/180)*a3)))

define(cea0, calc(cos((pi/180)*ea0)))
define(cea1, calc(cos((pi/180)*ea1)))
define(cea2, calc(cos((pi/180)*ea2)))
define(cea3, calc(cos((pi/180)*ea3)))

define(sea0, calc(sin((pi/180)*ea0)))
define(sea1, calc(sin((pi/180)*ea1)))
define(sea2, calc(sin((pi/180)*ea2)))
define(sea3, calc(sin((pi/180)*ea3)))

// Inner square x and y position
// x
define(x00, calc(s+x0))
define(x01, calc(-1.0*s+x0))
define(x02, calc(-1.0*s+x0))
define(x03, calc(s+x0))
// y
define(y00, calc(-1.0*s+y0))
define(y01, calc(-1.0*s+y0))
define(y02, calc(s+y0))
define(y03, calc(s+y0))

// Circle x and y positions
// x
define(x10, calc(r*ca0+x0))
define(x11, calc(r*ca1+x0))
define(x12, calc(r*ca2+x0))
define(x13, calc(r*ca3+x0))
// y
define(y10, calc(r*sa0+y0))
define(y11, calc(r*sa1+y0))
define(y12, calc(r*sa2+y0))
define(y13, calc(r*sa3+y0))

// Inner square x and y position middle curvatures
// x
define(ex00, calc(sc+x0))
define(ex01, x0)
define(ex02, calc(-1.0*sc+x0))
define(ex03, x0)
// y
define(ey00, y0)
define(ey01, calc(-1.0*sc+y0))
define(ey02, y0)
define(ey03, calc(sc+y0))

// Circle x and y positions middle curvatures
// x
define(ex10, calc(r*cea0+x0))
define(ex11, calc(r*cea1+x0))
define(ex12, calc(r*cea2+x0))
define(ex13, calc(r*cea3+x0))
// y
define(ey10, calc(r*sea0+y0))
define(ey11, calc(r*sea1+y0))
define(ey12, calc(r*sea2+y0))
define(ey13, calc(r*sea3+y0))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    vert(0, 0, Zb) vlabel(s0b)  //0
    vert(0, 1, Zb) vlabel(s1b)  //1
    vert(0, 2, Zb) vlabel(s2b)  //2
    vert(0, 3, Zb) vlabel(s3b)  //3
    
    vert(1, 0, Zb) vlabel(r0b)  //4
    vert(1, 1, Zb) vlabel(r1b)  //5
    vert(1, 2, Zb) vlabel(r2b)  //6
    vert(1, 3, Zb) vlabel(r3b)  //7

    vert(0, 0, Zm) vlabel(s0m)  //8
    vert(0, 1, Zm) vlabel(s1m)  //9
    vert(0, 2, Zm) vlabel(s2m)  //10
    vert(0, 3, Zm) vlabel(s3m)  //11
    
    vert(1, 0, Zm) vlabel(r0m)  //12
    vert(1, 1, Zm) vlabel(r1m)  //13
    vert(1, 2, Zm) vlabel(r2m)  //14
    vert(1, 3, Zm) vlabel(r3m)  //15
    
    vert(0, 0, Zt) vlabel(s0t)  //16
    vert(0, 1, Zt) vlabel(s1t)  //17
    vert(0, 2, Zt) vlabel(s2t)  //18
    vert(0, 3, Zt) vlabel(s3t)  //19
    
    vert(1, 0, Zt) vlabel(r0t)  //20
    vert(1, 1, Zt) vlabel(r1t)  //21
    vert(1, 2, Zt) vlabel(r2t)  //22
    vert(1, 3, Zt) vlabel(r3t)  //23
);

blocks
(
    //block0
    hex2D0(s1, s0, s3, s2)
    square
    (Ns0 Ns0 Nz0)
    simpleGrading (1 1 1)
    
    //block1
    hex2D0(s0, r0, r3, s3)
    innerCircle
    (Ni0 Ns0 Nz0)
    simpleGrading (stretchXY 1 1)
    
    //block2
    hex2D0(s3, r3, r2, s2)
    innerCircle
    (Ni0 Ns0 Nz0)
    simpleGrading (stretchXY 1 1)
    
    //block3
    hex2D0(s2, r2, r1, s1)
    innerCircle
    (Ni0 Ns0 Nz0)
    simpleGrading (stretchXY 1 1)
    
    //block4
    hex2D0(s1, r1, r0, s0)
    innerCircle
    (Ni0 Ns0 Nz0)
    simpleGrading (stretchXY 1 1)

    //block01
    hex2D1(s1, s0, s3, s2)
    square
    (Ns1 Ns1 Nz1)
    simpleGrading (1 1 stretchZ)
    
    //block11
    hex2D1(s0, r0, r3, s3)
    innerCircle
    (Ni1 Ns1 Nz1)
    simpleGrading (stretchXY 1 stretchZ)
    
    //block21
    hex2D1(s3, r3, r2, s2)
    innerCircle
    (Ni1 Ns1 Nz1)
    simpleGrading (stretchXY 1 stretchZ)
    
    //block31
    hex2D1(s2, r2, r1, s1)
    innerCircle
    (Ni1 Ns1 Nz1)
    simpleGrading (stretchXY 1 stretchZ)
    
    //block41
    hex2D1(s1, r1, r0, s0)
    innerCircle
    (Ni1 Ns1 Nz1)
    simpleGrading (stretchXY 1 stretchZ)
);

edges
(
    //Circle edges
    arc r3b r0b evert(1, 0, Zb)
    arc r0b r1b evert(1, 1, Zb)
    arc r1b r2b evert(1, 2, Zb)
    arc r2b r3b evert(1, 3, Zb)
    
    arc r3m r0m evert(1, 0, Zm)
    arc r0m r1m evert(1, 1, Zm)
    arc r1m r2m evert(1, 2, Zm)
    arc r2m r3m evert(1, 3, Zm)

    arc r3t r0t evert(1, 0, Zt)
    arc r0t r1t evert(1, 1, Zt)
    arc r1t r2t evert(1, 2, Zt)
    arc r2t r3t evert(1, 3, Zt)
    
    arc s3b s0b evert(0, 0, Zb)
    arc s0b s1b evert(0, 1, Zb)
    arc s1b s2b evert(0, 2, Zb)
    arc s2b s3b evert(0, 3, Zb)

    arc s3m s0m evert(0, 0, Zm)
    arc s0m s1m evert(0, 1, Zm)
    arc s1m s2m evert(0, 2, Zm)
    arc s2m s3m evert(0, 3, Zm)
    
    arc s3t s0t evert(0, 0, Zt)
    arc s0t s1t evert(0, 1, Zt)
    arc s1t s2t evert(0, 2, Zt)
    arc s2t s3t evert(0, 3, Zt)
);

patches
(
    patch inletOutlet
    (
        btQuad0(r0, r3)
        btQuad0(r1, r0)
        btQuad0(r2, r1)
        btQuad0(r3, r2)
        btQuad1(r0, r3)
        btQuad1(r1, r0)
        btQuad1(r2, r1)
        btQuad1(r3, r2)
    )
    
    wall Terrain
    (
        bottomQuad(s3, s0, s1, s2)
        bottomQuad(s3, r3, r0, s0)
        bottomQuad(s2, r2, r3, s3)
        bottomQuad(s1, r1, r2, s2)
        bottomQuad(s0, r0, r1, s1)
    )
    
    patch top
    (
        topQuad(s3, s0, s1, s2)
        topQuad(s3, r3, r0, s0)
        topQuad(s2, r2, r3, s3)
        topQuad(s1, r1, r2, s2)
        topQuad(s0, r0, r1, s1)
    )
);

mergePatchPairs
(
);

