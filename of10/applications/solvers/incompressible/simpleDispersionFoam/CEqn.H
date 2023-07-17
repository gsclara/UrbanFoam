   // Passive scalar transport
   volScalarField DC("DC", turbulence->nut()/Sct+turbulence->nu()/Sc);

   fvScalarMatrix CEqn
   (
       fvm::ddt(C)
     + fvm::div(phi, C)
     - fvm::laplacian(DC, C)
    ==
       fvModels.source(C)
   );

   CEqn.relax();
   fvConstraints.constrain(CEqn);
   CEqn.solve();
   fvConstraints.constrain(C);
