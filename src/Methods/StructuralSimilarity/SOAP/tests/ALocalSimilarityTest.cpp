//
// Created by Michael Heuer on 09.05.18.
//

#include <gtest/gtest.h>
#include <ParticlesVector.h>
#include <MolecularSpectrum.h>
#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "ExpansionSettings.h"
#include "Environment.h"
#include "TestMolecules.h"
#include "NeighborhoodExpander.h"

class ALocalSimilarityTest : public ::testing::Test {
public:

    MolecularGeometry molecule = TestMolecules::CO2::nuclei;
    double eps = std::numeric_limits<double>::epsilon()*1e3;
    void SetUp() override {};
};
// CAREFUL, ExpansionSettings made persist for following tests.

TEST_F(ALocalSimilarityTest , GenericNormalization) {
    ParticleKit::create(molecule);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::typeAgnostic;

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e2, e2), 1.0, eps);
};

TEST_F(ALocalSimilarityTest , SameEnvironmentsOnDifferentCenters) {
    ParticleKit::create(molecule);

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e2),1.0, eps);
};

TEST_F(ALocalSimilarityTest , Cross) {
    ParticleKit::create(molecule);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::typeAgnostic;

    Environment e0(molecule, molecule.atoms()[0].position()); // C
    Environment e1(molecule, molecule.atoms()[1].position()); // O1
    Environment e2(molecule, molecule.atoms()[2].position()); // O2

    auto val = LocalSimilarity::kernel(e0, e1);
    ASSERT_LT(val,1.0);
    ASSERT_GT(val,0.0);

    auto val2 = LocalSimilarity::kernel(e0, e2);
    ASSERT_LE(val2,1.0);
    ASSERT_GE(val2,0.0);

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1), LocalSimilarity::kernel(e0, e2), eps);
};


TEST_F(ALocalSimilarityTest, TypeSpecificNormalization) {
    ParticleKit::create(molecule);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e1),1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e2, e2),1.0, eps);
};

TEST_F(ALocalSimilarityTest, SameEnvironmentOnDifferentCentersGeneric) {
    ParticleKit::create(molecule);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::typeAgnostic;

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e2),1.0, eps);
};

TEST_F(ALocalSimilarityTest, SameEnvironmentOnDifferentCentersTypeSpecific) {
    ParticleKit::create(molecule);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    Environment e1(molecule,molecule.atoms()[1].position());
    Environment e2(molecule,molecule.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e1, e2),1.0, eps);
};

TEST_F(ALocalSimilarityTest, IsolatedSpecies) {
    ParticleKit::create(TestMolecules::CO2::isolatedNuclei);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto isolated = TestMolecules::CO2::isolatedNuclei;
    Environment e0(isolated, isolated.atoms()[0].position());
    Environment e1(isolated, isolated.atoms()[1].position());
    Environment e2(isolated, isolated.atoms()[2].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1),1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e2, e0),1.0, eps);
};

TEST_F(ALocalSimilarityTest, H2sameCenter) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(TestMolecules::H2::ElectronsInCores::normal);

    auto H2 = TestMolecules::H2::ElectronsInCores::normal;

    Environment a0(H2, H2.atoms()[0].position());
    Environment e0(H2, H2.electrons()[0].position());

    ASSERT_NEAR(LocalSimilarity::kernel(a0, e0), 1.0, eps);
}

TEST_F(ALocalSimilarityTest, H2sameEnvironment) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    ParticleKit::create(TestMolecules::H2::ElectronsInCores::normal);

    auto H2 = TestMolecules::H2::ElectronsInCores::normal;

    Environment a0(H2, H2.atoms()[0].position());
    Environment e0(H2, H2.electrons()[0].position());

    ASSERT_NEAR(LocalSimilarity::kernel(a0, e0), 1.0, eps);
}

TEST_F(ALocalSimilarityTest, twoOppositeElectrons) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto eaeb = TestMolecules::twoElectrons::oppositeSpin;
    ParticleKit::create(eaeb);

    Environment e0(eaeb, eaeb.electrons()[0].position());
    Environment e1(eaeb, eaeb.electrons()[1].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1),0.0,eps);
}

TEST_F(ALocalSimilarityTest, twoOppositeElectronsReversedOrder) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto ebea = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(ebea);

    Environment e0(ebea, ebea.electrons()[0].position());
    Environment e1(ebea, ebea.electrons()[1].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1),0.0,eps);
}

TEST_F(ALocalSimilarityTest, twoOppositeSpinElectronsComparision) {
    auto mol1 = TestMolecules::twoElectrons::oppositeSpin;
    auto mol2 = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(mol1);

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    Environment mol1e0(mol1, mol1.electrons()[0].position());
    Environment mol1e1(mol1, mol1.electrons()[1].position());

    Environment mol2e0(mol2, mol2.electrons()[0].position());
    Environment mol2e1(mol2, mol2.electrons()[1].position());

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol1e1), 0.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol2e1), 0.0, eps);

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol2e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol1e1), 1.0, eps);
};

TEST_F(ALocalSimilarityTest, twoOppositeSpinElectronsComparisionMs) {
    auto mol1 = TestMolecules::twoElectrons::oppositeSpin;
    auto mol2 = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(mol1);

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    MolecularSpectrum ms1(mol1);
    MolecularSpectrum ms2(mol2);

    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[NumberedType<int>(int(Spin::alpha), 0)],
                                        ms2.molecularCenters_[NumberedType<int>(int(Spin::alpha), 0)]), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[NumberedType<int>(int(Spin::beta), 0)],
                                        ms2.molecularCenters_[NumberedType<int>(int(Spin::beta), 0)]), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[NumberedType<int>(int(Spin::alpha), 0)],
                                        ms2.molecularCenters_[NumberedType<int>(int(Spin::beta), 0)]), 0.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(ms1.molecularCenters_[NumberedType<int>(int(Spin::beta), 0)],
                                        ms2.molecularCenters_[NumberedType<int>(int(Spin::alpha), 0)]), 0.0, eps);
}

TEST_F(ALocalSimilarityTest, twoAlphaElectrons) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto eaea = TestMolecules::twoElectrons::sameSpinAlpha;
    ParticleKit::create(eaea);

    Environment e0(eaea, eaea.electrons()[0].position());
    Environment e1(eaea, eaea.electrons()[1].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1),1.0,eps);
}


TEST_F(ALocalSimilarityTest, twoBetaElectrons) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto ebeb = TestMolecules::twoElectrons::sameSpinBeta;
    ParticleKit::create(ebeb);

    Environment e0(ebeb, ebeb.electrons()[0].position());
    Environment e1(ebeb, ebeb.electrons()[1].position());

    ASSERT_NEAR(LocalSimilarity::kernel(e0, e1),1.0,eps);
}

TEST_F(ALocalSimilarityTest, TypeSpecificAndAlchemicalComparison) {
    auto mol1 = TestMolecules::twoElectrons::oppositeSpin;
    auto mol2 = TestMolecules::twoElectrons::oppositeSpinReversedOrder;
    ParticleKit::create(mol1);

    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    Environment mol1e0(mol1, mol1.electrons()[0].position());
    Environment mol1e1(mol1, mol1.electrons()[1].position());

    Environment mol2e0(mol2, mol2.electrons()[0].position());
    Environment mol2e1(mol2, mol2.electrons()[1].position());

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol1e1), 0.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol2e1), 0.0, eps);

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol2e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol1e1), 1.0, eps);


    SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
    auto simMol1e0e1 = LocalSimilarity::kernel(mol1e0, mol1e1);
    auto simMol2e0e1 = LocalSimilarity::kernel(mol2e0, mol2e1);

    ASSERT_GT(simMol1e0e1, 0.0);
    ASSERT_LT(simMol1e0e1, 1.0);

    ASSERT_GT(simMol2e0e1, 0.0);
    ASSERT_LT(simMol2e0e1, 1.0);

    ASSERT_EQ(simMol1e0e1,simMol2e0e1);

    ASSERT_NEAR(LocalSimilarity::kernel(mol1e0, mol2e1), 1.0, eps);
    ASSERT_NEAR(LocalSimilarity::kernel(mol2e0, mol1e1), 1.0, eps);

};


TEST_F(ALocalSimilarityTest, DissociationIntoTwoIsolatedSpecies) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
    Cutoff::settings.radius = 2;// bohr
    Cutoff::settings.width = 1;// bohr // the inner plateau ends at 1
    ParticleKit::create({},{2,0}); //the particle kit consists of two alpha electrons

    // alter !both! environments by moving the second electron
    unsigned steps = 5  ;
    double rmax = 3; // bohr (cutoff is at 8 bohr)

    for (unsigned i = 0; i < steps; ++i) {
        double r = rmax*double(i)/double(steps-1);

        MolecularGeometry mol= {AtomsVector(),
                                ElectronsVector({{Spin::alpha, {0,0,0}},
                                                 {Spin::alpha, {0,0,r}}
                                                 })};

        Environment e0(mol, mol.electrons()[0].position());
        Environment e1(mol, mol.electrons()[1].position());

        ASSERT_EQ(LocalSimilarity::kernel(e0,e0),1);
        ASSERT_EQ(LocalSimilarity::kernel(e1,e1),1);
        ASSERT_EQ(LocalSimilarity::kernel(e0,e1),1);
    }
};

TEST_F(ALocalSimilarityTest, DissociationIntoOneIsolatedSpecies) {
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::typeAgnostic;
    Cutoff::settings.radius = 2;// bohr
    Cutoff::settings.width = 1;// bohr
    ParticleKit::create({},{3,0}); //the particle kit consists of three alpha electrons

    MolecularGeometry mol;
    mol= {AtomsVector(),ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,0}},{Spin::alpha, {0,0,-0.1}}})};
    Environment e0(mol, mol.electrons()[0].position());
    Environment e1(mol, mol.electrons()[1].position());
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e0),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e1,e1),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e1),1,eps);

    mol= {AtomsVector(),ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,0.5}},{Spin::alpha, {0,0,-0.1}}})};
    e0 = Environment(mol, mol.electrons()[0].position());
    e1 = Environment(mol, mol.electrons()[1].position());
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e0),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e1,e1),1,eps);
    ASSERT_GT(LocalSimilarity::kernel(e0,e1),0);
    ASSERT_LT(LocalSimilarity::kernel(e0,e1),1);

    mol= {AtomsVector(),ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,1.5}},{Spin::alpha, {0,0,-0.1}}})};
    e0 = Environment(mol, mol.electrons()[0].position());
    e1 = Environment(mol, mol.electrons()[1].position());
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e0),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e1,e1),1,eps);
    ASSERT_GT(LocalSimilarity::kernel(e0,e1),0);
    ASSERT_LT(LocalSimilarity::kernel(e0,e1),1);

    mol= {AtomsVector(),ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,2.0}},{Spin::alpha, {0,0,-0.1}}})};
    e0 = Environment(mol, mol.electrons()[0].position());
    e1 = Environment(mol, mol.electrons()[1].position());
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e0),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e1,e1),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e1),0,eps);

    mol= {AtomsVector(),ElectronsVector({{Spin::alpha,{0,0,0}},{Spin::alpha,{0,0,2.5}},{Spin::alpha, {0,0,-0.1}}})};
    e0 = Environment(mol, mol.electrons()[0].position());
    e1 = Environment(mol, mol.electrons()[1].position());
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e0),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e1,e1),1,eps);
    ASSERT_NEAR(LocalSimilarity::kernel(e0,e1),0,eps);
};