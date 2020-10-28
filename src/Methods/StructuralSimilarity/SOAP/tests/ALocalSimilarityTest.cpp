// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <ParticlesVector.h>
#include <MolecularSpectrum.h>
#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "SOAPSettings.h"
#include "Environment.h"
#include "TestMolecules.h"
#include "NeighborhoodExpander.h"

using namespace SOAP;

class ALocalSimilarityTest : public ::testing::Test {
public:

    MolecularGeometry molecule = TestMolecules::CO2::nuclei;
    double eps = std::numeric_limits<double>::epsilon() * 1e3;

    void SetUp() override {
        spdlog::set_level(spdlog::level::off);
    };
};

TEST_F(ALocalSimilarityTest, GenericNormalization) {
    ParticleKit::create(molecule);
    General::settings.mode = General::Mode::typeAgnostic;

    Environment e1(molecule, EnumeratedType<int>(int(Element::O), 0));
    Environment e2(molecule, EnumeratedType<int>(int(Element::O), 1));

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e2, e2), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, SameEnvironmentsOnDifferentCenters) {
    ParticleKit::create(molecule);

    Environment e1(molecule, EnumeratedType<int>(int(Element::O), 0));
    Environment e2(molecule, EnumeratedType<int>(int(Element::O), 1));

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e2), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, TwoDifferentParticlesOnSameCenter) {
    const auto mol = TestMolecules::threeElectrons::ionic;
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(mol);
    Radial::settings.nmax = 2;
    Angular::settings.lmax = 2;

    Environment e1(mol, EnumeratedType<int>(int(Spin::alpha), 1));
    Environment e2(mol, EnumeratedType<int>(int(Spin::beta), 0));

    ASSERT_LT(LocalSimilarity::kernel(e1, e2), 1.0 - eps);
};

TEST_F(ALocalSimilarityTest, Cross) {
    ParticleKit::create(molecule);
    General::settings.mode = General::Mode::typeAgnostic;

    Environment e0(molecule, EnumeratedType<int>(int(Element::C), 0));
    Environment e1(molecule, EnumeratedType<int>(int(Element::O), 0));
    Environment e2(molecule, EnumeratedType<int>(int(Element::O), 1));

    auto val = LocalSimilarity::kernel(e0, e1);
    ASSERT_LT(val, 1.0);
    ASSERT_GT(val, 0.0);

    auto val2 = LocalSimilarity::kernel(e0, e2);
    ASSERT_LE(val2, 1.0);
    ASSERT_GE(val2, 0.0);

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), LocalSimilarity::kernel(e0, e2), eps);
};


TEST_F(ALocalSimilarityTest, TypeSpecificNormalization) {
    ParticleKit::create(molecule);
    General::settings.mode = General::Mode::chemical;

    Environment e1(molecule, EnumeratedType<int>(int(Element::O), 0));
    Environment e2(molecule, EnumeratedType<int>(int(Element::O), 1));

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e2, e2), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, SameEnvironmentOnDifferentCentersGeneric) {
    ParticleKit::create(molecule);
    General::settings.mode = General::Mode::typeAgnostic;

    Environment e1(molecule, EnumeratedType<int>(int(Element::O), 0));
    Environment e2(molecule, EnumeratedType<int>(int(Element::O), 1));

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e2), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, SameEnvironmentOnDifferentCentersTypeSpecific) {
    ParticleKit::create(molecule);
    General::settings.mode = General::Mode::chemical;

    Environment e1(molecule, EnumeratedType<int>(int(Element::O), 0));
    Environment e2(molecule, EnumeratedType<int>(int(Element::O), 1));

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e2), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, IsolatedSpecies) {
    ParticleKit::create(TestMolecules::CO2::isolatedNuclei);
    General::settings.mode = General::Mode::chemical;

    auto isolated = TestMolecules::CO2::isolatedNuclei;
    Environment e0(isolated, molecule.findEnumeratedTypeByIndex(0));
    Environment e1(isolated, molecule.findEnumeratedTypeByIndex(1));
    Environment e2(isolated, molecule.findEnumeratedTypeByIndex(2));

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e2, e0), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, H2sameCenter) {
    General::settings.mode = General::Mode::chemical;
    ParticleKit::create(TestMolecules::H2::ElectronsInCores::normal);

    auto H2 = TestMolecules::H2::ElectronsInCores::normal;

    Environment a0(H2, EnumeratedType<int>(int(Element::H), 0));
    Environment e0(H2, EnumeratedType<int>(int(Spin::alpha), 0));

    ASSERT_LT(LocalSimilarity::kernel(a0, e0), 1.0 - eps);
}

TEST_F(ALocalSimilarityTest, twoOppositeElectrons) {
    General::settings.mode = General::Mode::chemical;

    auto eaeb = TestMolecules::twoElectrons::oppositeSpin;
    ParticleKit::create(eaeb);

    Environment e0(eaeb, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment e1(eaeb, EnumeratedType<int>(int(Spin::beta), 0));

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), 0.0, eps);
}

TEST_F(ALocalSimilarityTest, twoOppositeElectronsReversedOrder) {
    General::settings.mode = General::Mode::chemical;

    auto ebea = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(ebea);

    Environment e0(ebea, EnumeratedType<int>(int(Spin::beta), 0));
    Environment e1(ebea, EnumeratedType<int>(int(Spin::alpha), 0));

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), 0.0, eps);
}

TEST_F(ALocalSimilarityTest, twoOppositeSpinElectronsComparision) {
    auto mol1 = TestMolecules::twoElectrons::oppositeSpin;
    auto mol2 = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(mol1);

    General::settings.mode = General::Mode::chemical;

    Environment mol1e0(mol1, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment mol1e1(mol1, EnumeratedType<int>(int(Spin::beta), 0));

    Environment mol2e0(mol2, EnumeratedType<int>(int(Spin::beta), 0));
    Environment mol2e1(mol2, EnumeratedType<int>(int(Spin::alpha), 0));

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol1e1), 0.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol2e1), 0.0, eps);

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol2e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol1e1), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, twoOppositeSpinElectronsComparisionMs) {
    auto mol1 = TestMolecules::twoElectrons::oppositeSpin;
    auto mol2 = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(mol1);

    General::settings.mode = General::Mode::chemical;

    MolecularSpectrum ms1(mol1);
    MolecularSpectrum ms2(mol2);

    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[EnumeratedType<int>(int(Spin::alpha), 0)],
                                        ms2.molecularCenters_[EnumeratedType<int>(int(Spin::alpha), 0)]), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[EnumeratedType<int>(int(Spin::beta), 0)],
                                        ms2.molecularCenters_[EnumeratedType<int>(int(Spin::beta), 0)]), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[EnumeratedType<int>(int(Spin::alpha), 0)],
                                        ms2.molecularCenters_[EnumeratedType<int>(int(Spin::beta), 0)]), 0.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[EnumeratedType<int>(int(Spin::beta), 0)],
                                        ms2.molecularCenters_[EnumeratedType<int>(int(Spin::alpha), 0)]), 0.0, eps);
}

TEST_F(ALocalSimilarityTest, twoAlphaElectrons) {
    General::settings.mode = General::Mode::chemical;

    auto eaea = TestMolecules::twoElectrons::sameSpinAlpha;
    ParticleKit::create(eaea);

    Environment e0(eaea, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment e1(eaea, EnumeratedType<int>(int(Spin::alpha), 1));

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), 1.0, eps);
}


TEST_F(ALocalSimilarityTest, twoBetaElectrons) {
    General::settings.mode = General::Mode::chemical;

    auto ebeb = TestMolecules::twoElectrons::sameSpinBeta;
    ParticleKit::create(ebeb);

    Environment e0(ebeb, EnumeratedType<int>(int(Spin::beta), 0));
    Environment e1(ebeb, EnumeratedType<int>(int(Spin::beta), 1));

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), 1.0, eps);
}

TEST_F(ALocalSimilarityTest, ChemicalAndAlchemicalComparison) {
    auto mol1 = TestMolecules::twoElectrons::oppositeSpin;
    auto mol2 = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(mol1);

    General::settings.mode = General::Mode::chemical;
    Environment mol1e0(mol1, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment mol1e1(mol1, EnumeratedType<int>(int(Spin::beta), 0));

    Environment mol2e0(mol2, EnumeratedType<int>(int(Spin::beta), 0));
    Environment mol2e1(mol2, EnumeratedType<int>(int(Spin::alpha), 0));

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol1e1), 0.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol2e1), 0.0, eps);

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol2e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol1e1), 1.0, eps);


    General::settings.pairSimilarities[{int(Spin::alpha), int(Spin::beta)}] = 0.5;
    General::settings.mode = General::Mode::alchemical;
    auto simMol1e0e1 = LocalSimilarity::kernel(mol1e0, mol1e1);
    auto simMol2e0e1 = LocalSimilarity::kernel(mol2e0, mol2e1);

    ASSERT_GT(simMol1e0e1, 0.0);
    ASSERT_LT(simMol1e0e1, 1.0);

    ASSERT_GT(simMol2e0e1, 0.0);
    ASSERT_LT(simMol2e0e1, 1.0);

    ASSERT_EQ(simMol1e0e1, simMol2e0e1);

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol2e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol1e1), 1.0, eps);

};


TEST_F(ALocalSimilarityTest, DissociationIntoTwoIsolatedSpecies) {
    General::settings.mode = General::Mode::chemical;
    Cutoff::settings.radius = 2;// bohr
    Cutoff::settings.width = 1;// bohr // the inner plateau ends at 1
    ParticleKit::create({}, {2, 0}); //the particle kit consists of two alpha electrons

    // alter !both! environments by moving the second electron
    unsigned steps = 5;
    double rmax = 3; // bohr (cutoff is at 8 bohr)

    for (unsigned i = 0; i < steps; ++i) {
        double r = rmax * double(i) / double(steps - 1);

        MolecularGeometry mol = {AtomsVector(),
                                 ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                                  {Spin::alpha, {0, 0, r}}
                                                 })};

        Environment e0(mol, EnumeratedType<int>(int(Spin::alpha), 0));
        Environment e1(mol, EnumeratedType<int>(int(Spin::alpha), 1));

        ASSERT_EQ(LocalSimilarity::kernel(e0, e0), 1);
        ASSERT_EQ(LocalSimilarity::kernel(e1, e1), 1);
        ASSERT_EQ(LocalSimilarity::kernel(e0, e1), 1);
    }
};

TEST_F(ALocalSimilarityTest, DissociationIntoOneIsolatedSpecies) {
    General::settings.mode = General::Mode::typeAgnostic;
    Cutoff::settings.radius = 2;// bohr
    Cutoff::settings.width = 1;// bohr
    ParticleKit::create({}, {3, 0}); //the particle kit consists of three alpha electrons

    MolecularGeometry mol;
    mol = {AtomsVector(), ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                           {Spin::alpha, {0, 0, 0}},
                                           {Spin::alpha, {0, 0, -0.1}}})};
    Environment e0(mol, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment e1(mol, EnumeratedType<int>(int(Spin::alpha), 1));
    ASSERT_NEAR(LocalSimilarity::kernel(e0, e0), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e1, e1), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), 1, eps);

    mol = {AtomsVector(), ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                           {Spin::alpha, {0, 0, 0.5}},
                                           {Spin::alpha, {0, 0, -0.1}}})};
    Environment e2(mol, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment e3(mol, EnumeratedType<int>(int(Spin::alpha), 1));
    ASSERT_NEAR(LocalSimilarity::kernel(e2, e2), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e3, e3), 1, eps);
    ASSERT_GT(LocalSimilarity::kernel(e2, e3), 0);
    ASSERT_LT(LocalSimilarity::kernel(e2, e3), 1);

    mol = {AtomsVector(), ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                           {Spin::alpha, {0, 0, 1.5}},
                                           {Spin::alpha, {0, 0, -0.1}}})};
    Environment e4(mol, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment e5(mol, EnumeratedType<int>(int(Spin::alpha), 1));
    ASSERT_NEAR(LocalSimilarity::kernel(e4, e4), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e5, e5), 1, eps);
    ASSERT_GT(LocalSimilarity::kernel(e4, e5), 0);
    ASSERT_LT(LocalSimilarity::kernel(e4, e5), 1);

    mol = {AtomsVector(), ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                           {Spin::alpha, {0, 0, 2.0}},
                                           {Spin::alpha, {0, 0, -0.1}}})};
    Environment e6(mol, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment e7(mol, EnumeratedType<int>(int(Spin::alpha), 1));
    ASSERT_NEAR(LocalSimilarity::kernel(e6, e6), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e7, e7), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e6, e7), 0, eps);

    mol = {AtomsVector(), ElectronsVector({{Spin::alpha, {0, 0, 0}},
                                           {Spin::alpha, {0, 0, 2.5}},
                                           {Spin::alpha, {0, 0, -0.1}}})};
    Environment e8(mol, EnumeratedType<int>(int(Spin::alpha), 0));
    Environment e9(mol, EnumeratedType<int>(int(Spin::alpha), 1));
    ASSERT_NEAR(LocalSimilarity::kernel(e8, e8), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e9, e9), 1, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e8, e9), 0, eps);
};