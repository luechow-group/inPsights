//
// Created by heuer on 12.12.18.
//

#ifndef INPSIGHTS_GLOBALPERMUTATIONSORTER_H
#define INPSIGHTS_GLOBALPERMUTATIONSORTER_H

#include "Sample.h"
#include "SimilarReferences.h"
#include <StructuralSimilarity.h>
#include <spdlog/spdlog.h>

class GlobalPermutationSorter{
public:
    GlobalPermutationSorter(
            const AtomsVector& atoms,
            std::vector<Sample> &samples,
            std::vector<std::vector<SimilarReferences>> &globallyClusteredMaxima,
            std::vector<std::vector<SimilarReferences>> &globallyPermutationallyInvariantClusteredMaxima)
            :
            atoms_(atoms),
            samples_(samples),
            globallyClusteredMaxima_(globallyClusteredMaxima),
            globallyPermutationallyInvariantClusteredMaxima_(globallyPermutationallyInvariantClusteredMaxima)
    {
        ParticleKit::create(atoms_, (*samples.begin()).sample_);
    };

    void sort(){
        // compare all electron positions regardless of spin
        ExpansionSettings::defaults();
        //ExpansionSettings::Radial::nmax = 3;
        //ExpansionSettings::Angular::lmax = 3;
        ExpansionSettings::mode = ExpansionSettings::Mode::typeAgnostic;

        auto numberOfClusters = globallyClusteredMaxima_.size();

        std::vector<std::pair<std::vector<std::vector<SimilarReferences>>::iterator, MolecularSpectrum>>
        spectra(numberOfClusters);

        std::vector<std::pair<std::vector<std::vector<SimilarReferences>>::iterator, MolecularSpectrum>>
                identicalSpectra(0);

        spdlog::info("Spectra to calculate: {}", numberOfClusters);

        //collect representative structures before loop
        //std::vector<ElectronsVector> representativeClusterStructures(numberOfClusters);
        //for (size_t i = 0; i < globallyClusteredMaxima_.size(); ++i)
        //    representativeClusterStructures[i] = globallyClusteredMaxima_[i][0].representativeReference().maximum();

        double start = omp_get_wtime();
        //#pragma acc data copy(atomsVector,electronsVectors) create(spectra)
        //#pragma acc kernels

        #pragma omp parallel for default(none) shared(start, spectra, atoms_, globallyClusteredMaxima_)
        for (auto it = globallyClusteredMaxima_.begin(); it < globallyClusteredMaxima_.end(); ++it) {
            auto index = std::distance(globallyClusteredMaxima_.begin(),it);
            spectra[std::distance(globallyClusteredMaxima_.begin(),it)] =
                    {it, MolecularSpectrum({atoms_,(*it)[0].representativeReference().maximum()})};
            spdlog::info("Thread {0} wrote element i={1}, elapsed time: {2}",
                    omp_get_thread_num(), index, omp_get_wtime()-start);
        }

        assert(!spectra.empty() && "The vector cannot be empty");

        double identityThreshold = 0.98;

        globallyPermutationallyInvariantClusteredMaxima_.resize(0);


        //spdlog::info("Size before  = {}", globallyClusteredMaxima_.size());
        //spdlog::info("Size before  = {}", globallyPermutationallyInvariantClusteredMaxima_.size());
        //globallyPermutationallyInvariantClusteredMaxima_.emplace_back(std::move(*globallyClusteredMaxima_.begin()));
        //globallyClusteredMaxima_.erase(globallyClusteredMaxima_.begin());

        identicalSpectra.emplace_back(std::move(*spectra.begin()));
        spectra.erase(spectra.begin());

        //spdlog::info("Size after  = {}", globallyClusteredMaxima_.size());
        //spdlog::info("Size after  = {}", globallyPermutationallyInvariantClusteredMaxima_.size());



        auto it =  spectra.begin();
        while (it !=  spectra.end()) {

            // compare spectra
            //  first spectrum of

            // compare with all in list
            for(auto jt = identicalSpectra.begin(); jt != identicalSpectra.end(); jt++) {
                auto kdist = StructuralSimilarity::kernel((*it).second,(*jt).second);
                if (kdist <= identityThreshold) {

                    auto permutationallySimilarReferences = jt.base()->first.base();
                    // apply transformation that accounts for spatial symmetry to all similar refs and contained samples

                }
            }
        }


        //for (size_t i = 0; i < numberOfClusters; ++i) {
        //    //std::vector<std::vector<unsigned >>::iterator it;
        //    //printf("i=%d\n",i);
        //    for (size_t j = i+1; j < numberOfClusters; ++j) {
        //        auto kdist = StructuralSimilarity::kernel(spectra[i], spectra[j]);
        //        if(kdist >= threshold) {
        //            (*it).push_back(i);
        //            break;
        //        }
        //    }
        //    //// didn't fit into any other cluster
        //    //if(it == clusters.end()) {
        //    //    clusters.push_back({i});
        //    //}
        //}


    }

private:
    AtomsVector atoms_;
    std::vector<Sample> &samples_;
    std::vector<std::vector<SimilarReferences>> &globallyClusteredMaxima_;
    std::vector<std::vector<SimilarReferences>> &globallyPermutationallyInvariantClusteredMaxima_;
};

#endif //INPSIGHTS_GLOBALPERMUTATIONSORTER_H
