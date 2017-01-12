#include <iostream>

/*TODO Remember to set the amolqc path for runtime in the configuration file
 ** AMOLQC=/Users/michaelheuer/amolqcGUI/src/LibAmolqc/amolqc/
 *TODO Remember to put the t.wf file in the executable folders
 **    /Users/michaelheuer/amolqcGUI/build/cmake-build-debug/src/LibAmolqc/t.wf
 **    /Users/michaelheuer/amolqcGUI/build/cmake-build-release/src/LibAmolqc/t.wf
 */

typedef enum {GAUSSIAN, DENSITY, LMO} initPosType;

extern "C" {
   void amolqc_init();
   void amolqc_set_wf(int* nElecs, int* nAtoms);
   void amolqc_initial_positions(initPosType mode, int nElecs, double x[]);
   void amolqc_eloc(double x[], int n, double* phi, double* u, double grad[], double* elocal);
}

int main(int argc, char const *argv[])
{
   int nElecs, nAtoms;
   double Phi, U, E_local;


   amolqc_init();
   amolqc_set_wf(&nElecs, &nAtoms);
   std::cout << " # of elecs: " << nElecs << "  # of atoms: " << nAtoms << std::endl;

   double* x = new double[3*nElecs];
   double* drift = new double[3*nElecs];

   amolqc_initial_positions(DENSITY, nElecs, x);

   std::cout << " random electron positions (created with density mode): " << std::endl;
   for (int i = 0; i < nElecs; i++) {
      std::cout << i << " " << x[3*i] << " " << x[3*i+1] << " " << x[3*i+2] << std::endl;
   }

   amolqc_eloc(x, nElecs, &Phi, &U, drift, &E_local);
   std::cout << " Phi = " << Phi << "  U = " << U << " E_local = " << E_local << std::endl;
   std::cout << " drift = grad Psi / Psi = " << std::endl;

   for (int i = 0; i < nElecs; i++) {
      std::cout << i << " " << drift[3*i] << " " << drift[3*i+1] << " " << drift[3*i+2] << std::endl;
   }

   delete[] x;
   delete[] drift;

   return 0;
}

