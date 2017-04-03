/* minimal amolqc C++ interoperability test

   compile with GNU Compiler:
   g++ -c cpp_t.cpp
   link with g++ requires fortran runtime (-lgfortran) and openmp
   g++ -fopenmp -g -o cpp_t cpp_t.o src/libamolqc.a utils/libutils.a -Lutils -lutils 
      -L/usr/local/lib -llapack -lblas -lgfortran 
   link with gfortran requires C++ runtime lib:
   gfortran -fopenmp -g -o cpp_t cpp_t.o src/libamolqc.a utils/libutils.a -Lutils -lutils 
      -L/usr/local/lib -llapack -lblas -lstdc++   
   
   example assumes installation of compiler in /usr/local/bin, and libraries in /usr/local/lib
   do not use the preinstalled clang version of g++ on OSX

   Arne Luechow, 2017
 */

#include <iostream>

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

