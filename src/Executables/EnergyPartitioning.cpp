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



int main(int argc, char *argv[]) {

    Logger::initialize();
    auto console = spdlog::get(Logger::name);
    double identicalDistThresh = 0.01;
    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;
    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("Ethane.bin");
    console->info("number of refs {}",globallyIdenticalMaxima.size());
    auto numberOfElectrons = globallyIdenticalMaxima[0].maximum_.numberOfEntities();


    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identicalDistThresh);
    globalIdentiySorter.sort();
    console->info("finished id sort");
    console->flush();


    console->info("total elems {}",std::distance(globallyIdenticalMaxima.begin(),globallyIdenticalMaxima.end()));


    double similarDistThresh = 0.1;
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(globallyIdenticalMaxima, similarReferencesVector,similarDistThresh);

    globalSimilaritySorter.sort();

    console->info("total elems {}",similarReferencesVector.size());


    //std::vector<SimilarReference> clusteredReferences;

    /*Clustering*/
    // make data
    /*std::vector<Eigen::VectorXd> data;
    for(auto& simRefVector : similarReferencesVector){
        data.push_back((*simRefVector.representativeReferenceIterator).maximum_.positionsVector().asEigenVector());
    }

    // needs bestMatchDistance
    DensityBasedScan<double> dbscan(data);

    auto nClusters = dbscan.findClusters(0.20001, 5);
    auto result = dbscan.getLabels();*/


    /*Energy Partitioning*/
    Statistics::RunningStatistics<Eigen::VectorXd> EkinStats;
    Statistics::RunningStatistics<Eigen::MatrixXd> EpotStats;


    Eigen::VectorXd ekin;
    Eigen::MatrixXd epot;
    size_t sampleId, count, simCount, totalCount;

    totalCount = 0;

    int i=0;
    for(auto& simRefVector : similarReferencesVector){
        console->info("{}: contained similar refs {}",i,simRefVector.similarReferences_.size());
        i++;
        simCount = 0;

        EkinStats.reset();
        EpotStats.reset();

        // Representative reference
        count = 1+(*simRefVector.representativeReferenceIterator).associatedSampleIds_.size();
        simCount += count;
        sampleId = (*simRefVector.representativeReferenceIterator).id_;


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


    // Visuals
    std::string wavefunctionFilename = "Ethane-em.wf";
    auto wf = ElectronicWaveFunction::getInstance(wavefunctionFilename);
    auto atoms = wf.getAtomsVector();

    // Visualization
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    MoleculeWidget moleculeWidget;
    Qt3DCore::QEntity *root = moleculeWidget.createMoleculeWidget();

    AtomsVector3D(root, atoms);

    auto ev1 = (*similarReferencesVector.at(1).representativeReferenceIterator).maximum_;
    auto perm = similarReferencesVector.at(1).similarReferences_.at(0).perm_;
    auto ev2 = (*similarReferencesVector.at(1).similarReferences_.at(0).it_).maximum_;
    ev2.permute(perm);

    ElectronsVector3D(root, atoms, ev1, true);
    ElectronsVector3D(root, atoms, ev2, true);

    return app.exec();


//TODO CHECK PERMUTATION
// CHECK ADDITIONAL +1
//TODO DBSCAN

};
