//
// Created by Michael Heuer on 25.05.18.
//
#include <StructuralSimilarity.h>
#include <AmolqcFileImport/RefFileImporter.h>
#include <omp.h>
#include <sstream>

int main(int argc, char *argv[]) {

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

    RefFileImporter importer(argv[1]);//"Ethane-max.ref");
    auto atomsVector = importer.getAtomsVector();
    auto numberOfSuperstructures = importer.numberOfSuperstructures();

    // Settings
    ExpansionSettings::defaults();
    ExpansionSettings::Radial::nmax = nmax;
    ExpansionSettings::Angular::lmax = lmax;
    ExpansionSettings::mode = ExpansionSettings::Mode::Alchemical;
    ParticleKit::create(atomsVector, importer.getMaximaStructure(1,1));

    std::vector<MolecularSpectrum> spectra;
    spectra.reserve(numberOfSuperstructures);

    //#pragma omp parallel for shared(spectra)
    for (unsigned long i = 0; i < numberOfSuperstructures; ++i) {
        MolecularGeometry mol(atomsVector,importer.getMaximaStructure(i+1,1));
        //std::cout << mol;
        spectra.emplace(spectra.begin()+i,mol);
        std::cout << "done" << std::endl;
    }

    std::cout << "{\n";
    for (unsigned long i = 0; i < spectra.size(); ++i) {
        const auto& A = spectra[i];
        std::cout << "{";
        for (unsigned long j = i+1; j < spectra.size(); ++j) {
            const auto& B = spectra[j];
            std::cout <<"{"<< i << ","<< j << ","<< StructuralSimilarity::kernel(A,B)<< "},";
        }
        std::cout << "},\n";
    }
    std::cout << "}";
}