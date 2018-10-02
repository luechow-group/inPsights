//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <DensityBasedScan.h>
#include <CoulombPotential.h>
#include <Statistics.h>
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
    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("Ethane.bin");
    console->info("number of inital refs {}",globallyIdenticalMaxima.size());

    double identicalDistThresh = 0.01;
    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identicalDistThresh);
    globalIdentiySorter.sort();
    console->info("number of elements after identity sort {}",globallyIdenticalMaxima.size());


    //Clustering
    /*DensityBasedScan<double, Reference, wrapper2> dbscan(globallyIdenticalMaxima);
    auto nClusters = dbscan.findClusters(0.1*2+0.01, 1); // why multiplication by 2 is needed?
    auto labels = dbscan.getLabels();
    console->info("number of clusters {}",nClusters);

    for (int label : labels) {
        std::cout << label << std::endl;
    }*/

    double similarDistThresh = 0.2;
    std::vector<SimilarReferences> globallySimilarMaxima;
    GlobalSimilaritySorter globalSimilaritySorter(globallyIdenticalMaxima, globallySimilarMaxima,similarDistThresh);
    globalSimilaritySorter.sort();
    console->info("number of elements after similarity sort {}",globallySimilarMaxima.size());

    //Clustering
    DensityBasedScan<double, SimilarReferences, wrapper> dbscan(globallySimilarMaxima);
    auto nClusters = dbscan.findClusters(similarDistThresh*2+0.01, 1); // why multiplication by 2 is needed?
    auto labels = dbscan.getLabels();
    console->info("number of clusters {}",nClusters);


    //unify similar refs in the same cluster

    // for all similar refs


    std::vector<std::vector<SimilarReferences>> clusteredGloballySimilarMaxima(nClusters);

    for (int i = 0; i < nClusters; ++i) {
        auto labelCount = std::count(labels.begin(), labels.end(), i);

        int foundCount = 0;

        for (auto it = globallySimilarMaxima.begin(); it != globallySimilarMaxima.end(); ++it) {
            auto label = labels[std::distance(globallySimilarMaxima.begin(),it)];
            if (label == i) {
                clusteredGloballySimilarMaxima[i].emplace_back(std::move(*it));
            }
        }
    }

    for (auto it = clusteredGloballySimilarMaxima.begin(); it != clusteredGloballySimilarMaxima.end(); ++it) {
        int count = 0;
        for (auto & simRef : (*it)) count += simRef.similarReferences_.size()+1;
        console->info("cluster: {}, count: {}",std::distance(clusteredGloballySimilarMaxima.begin(),it),count);
    }

/*
    //Energy Partitioning
    Statistics::RunningStatistics<Eigen::VectorXd> EkinStats;
    Statistics::RunningStatistics<Eigen::MatrixXd> EpotStats;


    Eigen::VectorXd ekin;
    Eigen::MatrixXd epot;
    size_t sampleId, count, simCount, totalCount;

    totalCount = 0;

    int i=0;
    for(auto& simRefVector : globallySimilarMaxima){
        console->info("{}: contained similar refs {}",i,simRefVector.similarReferences_.size());
        i++;
        simCount = 0;

        EkinStats.reset();
        EpotStats.reset();

        // Representative reference
        count = 1+(*simRefVector.repRefIt_).associatedSampleIds_.size();
        simCount += count;
        sampleId = (*simRefVector.repRefIt_).id_;


        ekin = samples[sampleId].kineticEnergies_;
        epot = CoulombPotential::energies(samples[sampleId].sample_);

        console->info("rep ref count {}" ,count);
        EkinStats.add(ekin,unsigned(count));
        EpotStats.add(epot,unsigned(count));

        // Iterate over references being similar to the representative reference.
        for(const auto& simRef : simRefVector.similarReferences_){


            //TODO CONTINUE HERE
            globallyIdenticalMaxima[0];

            // Hier werden irgendwo maxima counts vergessen#
            count = 1+(*simRef.it_).associatedSampleIds_.size();
            simCount += count;
            sampleId = (*simRef.it_).id_;

            auto sampleCopy = samples[sampleId].sample_;
            sampleCopy.permute(simRef.perm_);

            ekin = simRef.perm_ * (samples[sampleId].kineticEnergies_);
            epot = CoulombPotential::energies(sampleCopy);

            EkinStats.add(ekin,unsigned(count));
            EpotStats.add(epot,unsigned(count));
        }
        console->info("sim ref count {}",simCount);
        std::cout <<"mean: ("<<EkinStats.getWeightedSum()<<")\n"<< EkinStats.mean().transpose() << std::endl;
        if(simCount >=2)
            std::cout <<"stdv:\n"<< EkinStats.standardDeviation().transpose() << std::endl<< std::endl;

        totalCount += simCount;
    }
    console->info("overall count {}",totalCount);
*/

    // Visuals
    /*std::string wavefunctionFilename = "Acetone-em.wf";
    auto wf = ElectronicWaveFunction::getInstance(wavefunctionFilename);
    auto atoms = wf.getAtomsVector();

    // Visualization
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();
    AtomsVector3D(root, atoms);

    auto ev1 = (*globallySimilarMaxima.at(1).repRefIt_).maximum_;
    auto perm = globallySimilarMaxima.at(1).similarReferences_.at(0).perm_;
    auto ev2 = (*globallySimilarMaxima.at(1).similarReferences_.at(0).it_).maximum_;
    ev2.permute(perm);
    ElectronsVector3D(root, atoms, ev1, true);
    ElectronsVector3D(root, atoms, ev2, true);

    return app.exec();*/


//TODO CHECK PERMUTATION
// CHECK ADDITIONAL +1
//TODO DBSCAN

};
