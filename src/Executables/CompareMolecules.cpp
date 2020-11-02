// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>
#include <iomanip>
#include <yaml-cpp/yaml.h>
#include <MolecularGeometry.h>

#include <ParticleKit.h>
#include <SOAPSettings.h>
#include <MolecularSpectrum.h>
#include <StructuralSimilarity.h>
#include <LocalSimilarity.h>

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.yml " << "clusterID" << std::endl;
        return 1;
    }
    std::string inputFilename = argv[1];

    YAML::Node doc = YAML::LoadFile(argv[1]);

    auto soapSettingsNode = doc["SOAP"];
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
    std::cout << node << std::endl;

    AtomKit akit;
    ElectronKit ekit;
    std::vector<MolecularGeometry> mols;

    std::vector<AtomsVector> atomsCollection;

    for(auto geomNode : doc["Geometries"]) {
        AtomsVector atoms;
        ElectronsVector electrons;

        if(geomNode["File"].IsDefined()) {
            auto molDoc = YAML::LoadAllFromFile(geomNode["File"]["name"].as<std::string>())[1];
            atoms = molDoc["Atoms"].as<AtomsVector>();
            auto ids = geomNode["File"]["ids"];
            
            for (auto id : ids) {
                auto clusterId = id["clusterId"].as<unsigned>();
                auto structureId = id["structureId"].as<unsigned>();
                electrons = molDoc["Clusters"][clusterId]["Structures"][structureId].as<ElectronsVector>();
                //electrons = molDoc["Clusters"][clusterId]["SampleAverage"].as<ElectronsVector>();
                mols.emplace_back(MolecularGeometry({}, electrons));
                //mols.emplace_back(MolecularGeometry(atoms, electrons));
                atomsCollection.emplace_back(atoms);
            }
        } else if (geomNode["Manual"].IsDefined()){
            atoms = geomNode["Manual"]["Atoms"].as<AtomsVector>();
            electrons = geomNode["Manual"]["Electrons"].as<ElectronsVector>();

            mols.emplace_back(MolecularGeometry({}, electrons));
            //mols.emplace_back(MolecularGeometry(atoms, electrons));
            atomsCollection.emplace_back(atoms);
        }

        akit = SOAP::ParticleKit::merge(
                akit,
                SOAP::ParticleKit::internal::createAtomKitFromAtomsVector(atoms));
        ekit = SOAP::ParticleKit::merge(
                ekit,
                SOAP::ParticleKit::internal::createElectronKitFromElectronsVector(electrons));
    }

    SOAP::ParticleKit::create(akit, ekit);

    std::cout << std::setprecision(2) << std::endl;
   for (std::vector<MolecularGeometry>::size_type i = 0; i < mols.size()-1; ++i) {
        for (std::vector<MolecularGeometry>::size_type j = i+1; j < mols.size(); ++j) {
            std::cout << i << " " << j << " "
                      << SOAP::LocalSimilarity::kernel(
                              SOAP::Environment(mols[i],(
                                      atomsCollection[i][1].position()+atomsCollection[i][0].position()) / 2.0),
                              SOAP::Environment(mols[j],(
                                      atomsCollection[j][0].position()+atomsCollection[j][1].position()) / 2.0)
                              ) << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "Glob" << std::endl;

    std::vector<SOAP::MolecularSpectrum> specs;
    for(const auto& mol : mols)
        specs.emplace_back(mol);

    for (std::vector<MolecularGeometry>::size_type i = 0; i < specs.size()-1; ++i) {
        for (std::vector<MolecularGeometry>::size_type j = i+1; j < specs.size(); ++j) {
            std::cout << i << " " << j << " "
            << SOAP::StructuralSimilarity::kernel(specs[i],specs[j]) << std::endl;
        }
        std::cout << std::endl;
    }
}