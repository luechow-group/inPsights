// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>
#include <yaml-cpp/yaml.h>
#include <MolecularGeometry.h>
#include <StructuralSimilarity.h>
#include <ParticleKit.h>
#include <SOAPSettings.h>

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

    auto geometriesNode = doc["Geometries"];

    AtomKit akit;
    ElectronKit ekit;

    std::vector<MolecularGeometry> mols;
    for (const auto& node : geometriesNode) {
        auto mol = MolecularGeometry(
                node["Atoms"].as<AtomsVector>(),
                node["Electrons"].as<ElectronsVector>());

        akit = SOAP::ParticleKit::merge(
                akit,
                SOAP::ParticleKit::internal::createAtomKitFromAtomsVector(mol.atoms()));
        ekit = SOAP::ParticleKit::merge(
                ekit,
                SOAP::ParticleKit::internal::createElectronKitFromElectronsVector(mol.electrons()));

        mols.emplace_back(mol);
    }

    SOAP::ParticleKit::create(akit, ekit);

    for (std::vector<MolecularGeometry>::size_type i = 0; i < mols.size()-1; ++i) {
        for (std::vector<MolecularGeometry>::size_type j = i+1; j < mols.size(); ++j) {
            std::cout << i << " " << j << " "
            << SOAP::StructuralSimilarity::kernel(mols[i],mols[j]) << std::endl;
        }
    }

}