//
// Created by Michael Heuer on 20.06.18.
//

#include "SimpleSorter.h"
#include "MolecularSpectrum.h"

std::vector<std::vector<unsigned >> SimpleSorter::sort(std::vector<MolecularSpectrum> spectra, double threshold){
    assert(!spectra.empty() && "The vector cannot be empty");

    std::vector<std::vector<unsigned >> clusters;
    clusters.push_back({0});
    unsigned n = spectra.size();
    for (unsigned i = 1; i < n; ++i) {//iterate over all structures
        std::vector<std::vector<unsigned >>::iterator it;

        for (it = clusters.begin(); it != clusters.end(); ++it){
            //auto kdist = StructuralSimilarity::kernelDistance(spectra[i],spectra[spectrumIdxOfTheFirstElement]);
            auto kdist = StructuralSimilarity::kernel(spectra[i],spectra[(*it)[0]]);
            if(kdist >= threshold) {
                (*it).push_back(i);
                break;
            }
        }
        // didn't fit into any other cluster
        if(it == clusters.end()) {
            clusters.push_back({i});
        }
    }

    for (auto& cluster : clusters){
        for(auto& element : cluster){
            std::cout << element << ",";
        }
        std::cout << std::endl;
    }
    return clusters;
};