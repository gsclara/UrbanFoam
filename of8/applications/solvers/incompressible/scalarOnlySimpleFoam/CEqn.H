   volScalarField DC("DC", turbulence->nu()/Sc+turbulence->nut()/Sct);

   fvScalarMatrix CEqn
   (
	fvm::ddt(C)
       +fvm::div(phi,C)
       -fvm::laplacian(DC,C)
      ==
        fvOptions(C)
   );

   CEqn.relax();
   fvOptions.constrain(CEqn);
   CEqn.solve();
   fvOptions.correct(C);
