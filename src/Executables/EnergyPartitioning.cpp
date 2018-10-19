#include <utility>

//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <DensityBasedScan.h>
#include <EnergyCalculator.h>

#include <MoleculeWidget.h>
#include <ElectronicWaveFunction.h>
#include <AtomsVector3D.h>
#include <ElectronsVector3D.h>
#include <QApplication>
#include <algorithm>


double wrapper(const SimilarReferences& s1, const SimilarReferences& s2) {
    //TODO add value to the metric
    return Metrics::bestMatchNorm<Eigen::Infinity,2>(
            s1.representativeReference().maximum(),
            s2.representativeReference().maximum());
};

struct SortElement{
    SortElement(
            std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch,
            std::vector<SimilarReferences>::iterator it)
            : bestMatch_(std::move(bestMatch)), it_(it)
    {}

    bool operator < (const SortElement& rhs) const {
        return bestMatch_.first < rhs.bestMatch_.first;
    }

    std::pair<double,Eigen::PermutationMatrix<Eigen::Dynamic>> bestMatch_;
    std::vector<SimilarReferences>::iterator it_;
};

int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    auto atoms = ElectronicWaveFunction::getInstance(std::string("Acetone-em.wf")).getAtomsVector();
    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("Acetone.bin");


    console->info("number of inital refs {}",globallyIdenticalMaxima.size());

    EnergyCalculator energyCalculator(samples,atoms);
    auto totalEnergies = energyCalculator.calculateTotalEnergies();

    console->info("Te= {}, Vee = {}, Ven = {}, Vnn = {}, Eges = {}",
            totalEnergies.Te,
            totalEnergies.Vee,
            totalEnergies.Ven,
            totalEnergies.Vnn,
            totalEnergies.totalEnergy());


    double identicalDistThresh = 0.01;
    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identicalDistThresh);
    globalIdentiySorter.sort();
    console->info("number of elements after identity sort {}",globallyIdenticalMaxima.size());

    double similarDistThresh = 0.2;
    std::vector<SimilarReferences> globallySimilarMaxima;
    GlobalSimilaritySorter globalSimilaritySorter(samples,globallyIdenticalMaxima, globallySimilarMaxima,similarDistThresh);
    globalSimilaritySorter.sort();
    console->info("number of elements after similarity sort {}",globallySimilarMaxima.size());

    //Clustering
    DensityBasedScan<double, SimilarReferences, wrapper> dbscan(globallySimilarMaxima);
    auto nClusters = dbscan.findClusters(similarDistThresh*2+0.01, 1); // why multiplication by 2 is needed?
    //TODO Permutations must be stored! Own DBSCAN implementation?
    auto labels = dbscan.getLabels();
    console->info("number of clusters {}",nClusters);

    std::vector<std::vector<SimilarReferences>> clusteredGloballySimilarMaxima(nClusters);

    for (int i = 0; i < nClusters; ++i) {
        for (auto it = globallySimilarMaxima.begin(); it != globallySimilarMaxima.end(); ++it) {
            auto label = labels[std::distance(globallySimilarMaxima.begin(),it)];
            if (label == i) {
                clusteredGloballySimilarMaxima[i].emplace_back(std::move(*it));
            }
        }
    }

    // orderd clusters
    for (auto& cluster : clusteredGloballySimilarMaxima) {
        std::sort(cluster.begin(), cluster.end());

        for (auto i = cluster.begin(); i != cluster.end(); ++i) {
            std::vector<SortElement> bestMatchDistances;

            for (auto j = i + 1; j != cluster.end(); ++j) {
                bestMatchDistances.emplace_back(
                        SortElement(Metrics::bestMatch<Eigen::Infinity, 2>(
                                j.base()->representativeReference().maximum(),
                                i.base()->representativeReference().maximum()), j)
                );
            }

            if (bestMatchDistances.size() > 1) {
                // permute and swap
                auto minIt = std::min_element(bestMatchDistances.begin(), bestMatchDistances.end());
                if (i + 1 != minIt.base()->it_)
                    std::iter_swap(i + 1, minIt.base()->it_);
            }
        }
    }

    //Statistics
    energyCalculator.calculateStatistics(clusteredGloballySimilarMaxima);
    std::ofstream yamlFile("energies.yml");
    yamlFile << energyCalculator.getYamlDocumentString();
    yamlFile.close();

    // Visuals
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();
    AtomsVector3D(root, atoms);

    for (auto i : clusteredGloballySimilarMaxima[35]){
        ElectronsVector3D(root, atoms, i.representativeReference().maximum(), false);
    }
    return QApplication::exec();

};
