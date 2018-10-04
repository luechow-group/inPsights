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
    return Metrics::bestMatchNorm<Eigen::Infinity,2>((*s1.repRefIt_).maximum_, (*s2.repRefIt_).maximum_);
};

int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    auto atoms = ElectronicWaveFunction::getInstance(std::string("Ethane-em.wf")).getAtomsVector();
    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("Ethane.bin");


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
    GlobalSimilaritySorter globalSimilaritySorter(globallyIdenticalMaxima, globallySimilarMaxima,similarDistThresh);
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

    //Statistics
    energyCalculator.calculateStatistics(clusteredGloballySimilarMaxima);


    // Visuals
    /*QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();
    AtomsVector3D(root, atoms);

    auto ev1 = (*clusteredGloballySimilarMaxima[0].at(0).repRefIt_).maximum_;
    auto ev2 = samples[(*clusteredGloballySimilarMaxima[0].at(0).repRefIt_).associatedSampleIds_[0]].sample_;

    //auto perm = globallySimilarMaxima.at(1).similarReferences_.at(0).perm_;
    //auto ev2 = (*clusteredGloballySimilarMaxima[0].at(0).similarReferences_.at(0).it_).maximum_;
    //ev2.permute(perm);
    ElectronsVector3D(root, atoms, ev1, true);
    ElectronsVector3D(root, atoms, ev2, true);

    return app.exec();*/

};
