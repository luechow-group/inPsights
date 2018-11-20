//
// Created by Michael Heuer on 25.05.18.
//
#include <omp.h>
#include <RefFileImporter.h>
#include <StructuralSimilarity.h>
#include <SimpleSorter.h>
#include <Serialization.h>

int main(int argc, char *argv[]) {

    const int maxNumThreads = omp_get_max_threads();
    printf("Threads used: %i\n", maxNumThreads);
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

    std::istringstream filenameInput(argv[1]);
    std::string filename;
    if (!(filenameInput >> filename)) {
        std::cerr << "Invalid filename " << argv[1] << '\n';
        return false;
    }

    std::istringstream nmaxInput(argv[2]);
    unsigned nmax;
    if (!(nmaxInput >> nmax) || nmax <= 0) {
        std::cerr << "Invalid number " << argv[2] << '\n';
        return false;
    }

    std::istringstream lmaxInput(argv[3]);
    unsigned lmax;
    if (!(lmaxInput >> lmax)){
        std::cerr << "Invalid number " << argv[3] << '\n';
        return false;
    }

    std::cout << "Filename: " << filename << std::endl;
    RefFileImporter importer(filename);
    auto atomsVector = importer.getAtomsVector();
    std::cout << atomsVector<< std::endl;

    // Settings
    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionSettings::Mode::alchemical;
    ExpansionSettings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 0.5;
    ParticleKit::create(atomsVector,importer.getMaximaStructure(1,1));
    std::cout << ExpansionSettings::toString() << "\n" << ParticleKit::toString() << std::endl;

    // Benchmark
    /*for (unsigned i = 1; i <= 14; ++i) {
        ExpansionSettings::Radial::nmax = i;
        ExpansionSettings::Angular::lmax = i;
        double start = omp_get_wtime();
        MolecularSpectrum ms1({atomsVector,importer.getMaximaStructure(1,1)});
        MolecularSpectrum ms2({atomsVector,importer.getMaximaStructure(2,1)});
        auto t1 = omp_get_wtime()-start;

        //auto k=LocalSimilarity::kernel({ms1.molecularCenters_[{6,0}]},{ms1.molecularCenters_[{6,1}]},ExpansionSettings::zeta);

        auto k = StructuralSimilarity::kernel(ms1,ms2);
        auto t2 = omp_get_wtime()-start;
        printf("{%d,%d,%.17g,%f,%f},\n",i,i, k, t1, t2);
    }*/

    unsigned numberOfSuperstructures = importer.numberOfSuperstructures();
    std::cout <<"Number of superstrcuctures in "<< filename << ": " << numberOfSuperstructures << std::endl;
    ElectronsVectorCollection electronsVectors;
    for (unsigned k = 0; k < numberOfSuperstructures; ++k) {
        electronsVectors.append(importer.getMaximaStructure(k+1,1));
    }
    std::cout << electronsVectors << std::endl;

    // Settings
    ExpansionSettings::defaults();
    ExpansionSettings::Radial::nmax = nmax;
    ExpansionSettings::Angular::lmax = lmax;
    ExpansionSettings::mode = ExpansionSettings::Mode::alchemical;
    ExpansionSettings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 0.5;
    ParticleKit::create(atomsVector,importer.getMaximaStructure(1,1));
    std::cout << ExpansionSettings::toString() << "\n" << ParticleKit::toString() << std::endl;

    std::vector<MolecularSpectrum> spectra(numberOfSuperstructures);
    printf("Spectra to calculate: %lu\n",spectra.size());

    double start = omp_get_wtime();
    //#pragma acc data copy(atomsVector,electronsVectors) create(spectra)
    //#pragma acc kernels
    #pragma omp parallel for default(none) shared(start,numberOfSuperstructures,spectra,atomsVector,electronsVectors)
    for (unsigned i = 0; i < numberOfSuperstructures; ++i) {
        spectra[i] = MolecularSpectrum({atomsVector,electronsVectors[i]});
        printf("Thread %d wrote element i=%ld\telapsed time: %fs\n", omp_get_thread_num(),i,omp_get_wtime()-start);
    }

    double simBorder=0.98;

    SimpleSorter sorter;
    auto clusters = sorter.sort(spectra,simBorder);
    std::cout << Serialization::yamlStringFrom(clusters,YAML::Flow) << std::endl;
    printf("elapsed time: %fs",omp_get_wtime()-start);
    // Inefficient
    /*
    Eigen::MatrixXd pairs = Eigen::MatrixXd::Identity(numberOfSuperstructures,numberOfSuperstructures);
    for (unsigned long i = 0; i < spectra.size(); ++i) {
        //#pragma omp parallel for default(none) shared(i,spectra,vec) firstprivate(ExpansionSettings::gamma)
        //#pragma acc kernels
        for (unsigned long j = i+1; j < spectra.size(); ++j) {
            pairs(i,j) = StructuralSimilarity::kernel(spectra[i],spectra[j],ExpansionSettings::gamma);
            printf("{%lu,%lu,%f} \n",i, j, pairs(i,j));
        }
    }
    printf("elapsed time: %fs",omp_get_wtime()-start);
*/




}
