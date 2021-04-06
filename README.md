# UrbanFoam
The repository contains different libraries and applications that can be used in urban simulations. OpenFOAM versions are 7 and 8, in their respective folders.

**OpenFOAM v7:**

A list of the different applications is: 
- (in development) applications/solvers/incompressible/simplez0Foam: it is a solver that attempts to take the z0 directly as a surface field based on architecture sermantics from 3D city models
- applications/utilities/surface/surfaceInfo: it is an application that simplifies surfaceCheck, but just printing the extreme coordinates of an .stl or .obj file without checking the surface properties or saving any additional .obj files
- applications/solvers/incompressible/simpleTFoam: it is a solver that includes the temperature field as a passive scalar 
- applications/solvers/heatTransfer/buoyantBoussinesqPimpleFoam: ported the solver from OF6
- applications/solvers/heatTransfer/buoyantBoussinesqSimpleFoam: ported the solver from OF6

A list of the different libraries is: 
- src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonz0WallFunction: wall function adapted from Parente et al. 2011b for atmospheric flows
- src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/nutz0WallFunction: wall function adapted from Parente et al. 2011b for atmospheric flows
- src/TurbulenceModels/incompressible/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatkz0WallFunction: ported the wall function from of-v2006 by ESI group, adapts the alphat for neutral ATMBL. NOTE: uses scalar instead of scalarField for now
- (in development) src/finiteVolume/fields/fvPatchFields/derived/timeVaryingInletOutletTke/timeVaryingInletOutletTke: time varying boundary condition for neutral atmospheric boundary layer turbulence kinetic energy with changing wind direction
- (in development) src/finiteVolume/fields/fvPatchFields/derived/timeVaryingInletOutletEpsilon/timeVaryingInletOutletEpsilon: time varying boundary condition for neutral atmospheric boundary layer epsilon with changing wind direction
- (in development) src/finiteVolume/fields/fvPatchFields/derived/timeVaryingInletOutletVelocity/timeVaryingInletOutletVelocity: time varying boundary condition for neutral atmospheric boundary layer velocity with changing wind direction

A list of the different functionObjects is:
- /src/functionObjects/fields/firstCellHeight: function object that calculates the height of the first cell next to a wall

**OpenFOAM v8:**

A list of the different applications is: 
- applications/solvers/heatTransfer/buoyantBoussinesqSimpleFoam: ported the solver from OF6; requires alphatJayatillekeWallFunction to be compiled

A list of the different libraries is: 
- src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatJayatillekeWallFunction: ported the wall function from OF6, it is used for incompressible heat transfer (buoyant) flows
- src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonz0WallFunction: wall function adapted from Parente et al. 2011b for atmospheric flows
- src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/nutz0WallFunction: wall function adapted from Parente et al. 2011b for atmospheric flows

A list of the different functionObjects is:
- /src/functionObjects/fields/firstCellHeight: function object that calculates the height of the first cell next to a wall
