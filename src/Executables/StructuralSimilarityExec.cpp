//
// Created by Michael Heuer on 25.05.18.
//
#include <StructuralSimilarity.h>
#include <AmolqcFileImport/RefFileImporter.h>
#include <omp.h>
#include <sstream>

int main(int argc, char *argv[]) {
    const int maxNumThreads = omp_get_max_threads();
    printf("Maximum number of threads for this machine: %i\n", maxNumThreads);
    omp_set_num_threads(maxNumThreads);

    if (argc < 4) {
        std::cout << "Arguments missing. \n"
                     "reference file argument must be a string ending with .ref\n"
                     "nmax argument must be a positive integer > 0\n"
                     "lmax argument must be a positive integer"
                  << std::endl;
        std::cout << "e.g. ./StructuralSimilarity.exe max.ref 10 10" << std::endl;
        return false;
    }

    std::istringstream nmaxString(argv[2]);
    unsigned nmax;
    if (!(nmaxString >> nmax) || nmax <= 0) {
        std::cerr << "Invalid number " << argv[2] << '\n';
        return false;
    }

    std::istringstream lmaxString(argv[3]);
    unsigned lmax;
    if (!(lmaxString >> lmax)){
        std::cerr << "Invalid number " << argv[3] << '\n';
        return false;
    }

    std::cout << argv[1] << std::endl;
    std::cout << "nmax: " << nmax << std::endl;
    std::cout << "lmax: " << lmax << std::endl;

    RefFileImporter importer(argv[1]);
    auto atomsVector = importer.getAtomsVector();
    auto numberOfSuperstructures = importer.numberOfSuperstructures();

    // Settings
    ExpansionSettings::defaults();
    ExpansionSettings::Radial::nmax = nmax;
    ExpansionSettings::Angular::lmax = lmax;
    ExpansionSettings::mode = ExpansionSettings::Mode::alchemical;
    ExpansionSettings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 1.0;
    ParticleKit::create(atomsVector, importer.getMaximaStructure(1,1));

    std::cout << ExpansionSettings::toString() << "\n" << ParticleKit::toString() << std::endl;

    std::vector<MolecularSpectrum> spectra(numberOfSuperstructures);
    printf("size at the beginning %lu\n",spectra.size());

    double start = omp_get_wtime();
    #pragma omp parallel for default(none) shared(start,numberOfSuperstructures,spectra,atomsVector,importer)
    for (unsigned long i = 0; i < numberOfSuperstructures; ++i) {
        spectra[i] = MolecularSpectrum({atomsVector,importer.getMaximaStructure(i+1,1)});
        printf("Thread %d wrote element i=%d\nelapsed time: %fs\n", omp_get_thread_num(),i,omp_get_wtime()-start);
    }
    //TODO RANDOMLY OCCURING ERROR
    //Assertion failed: (m <= substructuresData_[k].numberOfSubstructures_ && "m value must be smaller than mmax"), function calculateLine, file /Users/michaelheuer/Projects/Amolqcpp/src/AmolqcInterface/source/AmolqcFileImport/RefFileImporter.cpp, line 45.


    for (unsigned long i = 0; i < spectra.size(); ++i) {
        #pragma omp parallel for default(none) shared(i,spectra) firstprivate(ExpansionSettings::gamma)
        for (unsigned long j = i+1; j < spectra.size(); ++j) {
            double val = StructuralSimilarity::kernel(spectra[i],spectra[j],ExpansionSettings::gamma);
            printf("{%lu,%lu,%f} \n",i, j, val);
        }
    }
    printf("elapsed time: %fs",omp_get_wtime()-start);
}