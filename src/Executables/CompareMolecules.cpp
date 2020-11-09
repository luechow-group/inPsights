// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>
#include <iomanip>
#include <yaml-cpp/yaml.h>

#include <fstream>

#include <MolecularGeometry.h>
#include <Enumerate.h>
#include <ParticleKit.h>
#include <SOAPSettings.h>
#include <MolecularSpectrum.h>
#include <StructuralSimilarity.h>


struct LabeledMolecule {
    std::string label;
    MolecularGeometry molecule;
};

void writeResults(const std::string &inputFilename, const YAML::Emitter &emitter);

std::vector<LabeledMolecule> readMoleculesFromYaml(const YAML::Node &doc);

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.yml " << "clusterID" << std::endl;
        return 1;
    }
    std::string inputFilename = argv[1];

    YAML::Node inputYaml = YAML::LoadFile(argv[1]);

    auto soapSettingsNode = inputYaml["SOAP"];
    if (soapSettingsNode[SOAP::General::settings.name()])
        SOAP::General::settings = Settings::SOAP::General(soapSettingsNode);
    if (soapSettingsNode[SOAP::Radial::settings.name()])
        SOAP::Radial::settings = Settings::SOAP::Radial(soapSettingsNode);
    if (soapSettingsNode[SOAP::Angular::settings.name()])
        SOAP::Angular::settings = Settings::SOAP::Angular(soapSettingsNode);
    if (soapSettingsNode[SOAP::Cutoff::settings.name()])
        SOAP::Cutoff::settings = Settings::SOAP::Cutoff(soapSettingsNode);

    YAML::Node node;
    SOAP::General::settings.appendToNode(node);
    std::cout << inputYaml << std::endl;


    auto mols = readMoleculesFromYaml(inputYaml);

    std::vector<SOAP::MolecularSpectrum> specs;

    for(const auto& [i, mol] : enumerate(mols)) {
        std::cout << i+1 << " of " << mols.size() << std::endl;
        specs.emplace_back(mol.molecule);
    }

    using namespace YAML;
    Emitter emitter;
    emitter << BeginDoc << Comment("Used input from " + inputFilename + ".") << inputYaml << EndDoc;
    emitter << BeginDoc << BeginMap;

    std::cout << std::endl;
    for (std::vector<MolecularGeometry>::size_type i = 0; i < specs.size()-1; ++i) {
        emitter << Key << mols[i].label << Value <<  BeginMap;
        for (std::vector<MolecularGeometry>::size_type j = i+1; j < specs.size(); ++j) {
            auto value = SOAP::StructuralSimilarity::kernel(specs[i],specs[j]);
            std::cout << i << " " << j << " " << value << std::endl;
            emitter << Key << mols[j].label  << Value << value;
        }
        emitter << EndMap;
    }
    emitter << EndMap << EndDoc;


    writeResults(inputFilename, emitter);
}

std::vector<LabeledMolecule> readMoleculesFromYaml(const YAML::Node &doc) {
    std::vector<LabeledMolecule> mols;

    AtomKit akit;
    ElectronKit ekit;

    unsigned structureCounter = 0;
    for(auto geomNode : doc["Geometries"]) {
        AtomsVector atoms;
        ElectronsVector electrons;

        if(geomNode["File"].IsDefined()) {
            auto molDoc = YAML::LoadAllFromFile(geomNode["File"]["name"].as<std::string>())[1];
            atoms = {};//molDoc["Atoms"].as<AtomsVector>();
            auto ids = geomNode["File"]["ids"];

            for (auto id : ids) {
                auto clusterId = id["clusterId"].as<unsigned>();
                auto structureId = id["structureId"].as<unsigned>();

                std::string label = std::to_string(structureCounter);//std::to_string(clusterId) + "-" + std::to_string(structureId);
                if(id["label"].IsDefined())
                    label = id["label"].as<std::string>();

                electrons = molDoc["Clusters"][clusterId]["Structures"][structureId].as<ElectronsVector>();
                //electrons = molDoc["Clusters"][clusterId]["SampleAverage"].as<ElectronsVector>();

                mols.emplace_back<LabeledMolecule>({label,MolecularGeometry(atoms, electrons)});
                structureCounter++;
            }
        } else if (geomNode["Manual"].IsDefined()){
            atoms = geomNode["Manual"]["Atoms"].as<AtomsVector>();
            electrons = geomNode["Manual"]["Electrons"].as<ElectronsVector>();
            std::string label = std::to_string(structureCounter);

            if(geomNode["label"].IsDefined())
                label = geomNode["Manual"]["label"].as<std::string>();

            mols.emplace_back<LabeledMolecule>({label,MolecularGeometry(atoms, electrons)});
            structureCounter++;
        }

        akit = SOAP::ParticleKit::merge(
                akit,
                SOAP::ParticleKit::internal::createAtomKitFromAtomsVector(atoms));
        ekit = SOAP::ParticleKit::merge(
                ekit,
                SOAP::ParticleKit::internal::createElectronKitFromElectronsVector(electrons));
    }

    SOAP::ParticleKit::create(akit, ekit);

    return mols;
}

void writeResults(const std::string &inputFilename, const YAML::Emitter &emitter) {
    std::string resultsBaseName = inputFilename.substr(0, inputFilename.find('.')) + "-out";
    std::string resultsFilename = resultsBaseName + ".yml";
    if (std::ifstream(resultsFilename).good()){
        for (int i = 1; i < 100; i++){
            resultsFilename = resultsBaseName + "-" + std::to_string(i) + ".yml";
            if (not std::ifstream(resultsFilename).good()){
                break;
            }
        }
    }
    spdlog::info("Writing results into file \"{}\"", resultsFilename);

    auto outstring = emitter.c_str();
    spdlog::info("\n{}",outstring);

    std::ofstream yamlFile(resultsFilename);
    yamlFile << outstring;
    yamlFile.close();

    std::cout << "Wrote results to " << resultsFilename << std::endl;
}
