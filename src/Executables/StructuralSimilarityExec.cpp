//
// Created by Michael Heuer on 25.05.18.
//
#include <StructuralSimilarity.h>
#include <AmolqcFileImport/RefFileImporter.h>

int main(int argc, char *argv[]) {
    RefFileImporter importer("Ethane-max.ref");
    auto atomsVector = importer.getAtomsVector();
    auto numberOfSuperstructures = importer.numberOfSuperstructures();

    // Settings
    ExpansionSettings::defaults();
    ExpansionSettings::Radial::nmax = 5;
    ExpansionSettings::Angular::lmax = 5;
    ExpansionSettings::mode = ExpansionSettings::Mode::Alchemical;
    ParticleKit::create(atomsVector, importer.getMaximaStructure(1,1));

    std::vector<MolecularSpectrum> spectra;
    spectra.reserve(numberOfSuperstructures);

    for (unsigned long i = 0; i < numberOfSuperstructures; ++i) {
        MolecularGeometry mol(atomsVector,importer.getMaximaStructure(i+1,1));
        std::cout << mol;
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