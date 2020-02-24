# UrbanFoam
The repository contains different libraries and applications that can be used in urban simulations. 
A list of the different applications is: 
- (in development) applications/solvers/incompressible/simplez0Foam: it is a solver that attempts to take the z0 directly as a surface field based on architecture sermantics from 3D city models
- applications/utilities/surface/surfaceInfo: it is an application that simplifies surfaceCheck, but just printing the extreme coordinates of an .stl or .obj file without checking the surface properties or saving any additional .obj files

A list of the different libraries is: 
- src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/epsilonz0WallFunction: wall function adapted from Parente et al. 2011b for atmospheric flows
- src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/epsilonWallFunctions/nutz0WallFunction: wall function adapted from Parente et al. 2011b for atmospheric flows
- src/finiteVolume/fields/fvPatchFields/derived/timeVaryingInletOutletTke/timeVaryingInletOutletTke: time varying boundary condition for neutral atmospheric boundary layer turbulence kinetic energy with changing wind direction
- src/finiteVolume/fields/fvPatchFields/derived/timeVaryingInletOutletEpsilon/timeVaryingInletOutletEpsilon: time varying boundary condition for neutral atmospheric boundary layer epsilon with changing wind direction
- src/finiteVolume/fields/fvPatchFields/derived/timeVaryingInletOutletVelocity/timeVaryingInletOutletVelocity: time varying boundary condition for neutral atmospheric boundary layer velocity with changing wind direction
