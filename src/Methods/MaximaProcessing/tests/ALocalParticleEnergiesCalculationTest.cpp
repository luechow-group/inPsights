/* Copyright (C) 2020 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <gmock/gmock.h>
#include <LocalParticleEnergiesCalculation.h>
#include <CoulombPotential.h>
#include <Sample.h>
#include <Reference.h>
#include <Group.h>
#include <yaml-cpp/yaml.h>

using namespace testing;

class ALocalParticleEnergiesCalculationTest : public ::testing::Test {
public:

    void SetUp() override {
        //spdlog::set_level(spdlog::level::debug);
    }


    double sumAll(const LocalParticleEnergiesCalculator::LocalBondEnergyResults &energies) const {
        double sum = (energies.intraBond.mean() + energies.intraRest.mean() + energies.interBondRest.mean())(0, 0);

        for (Eigen::Index k = 0; k < energies.intraCores.mean().rows(); ++k)
            sum += energies.intraCores.mean()[k];

        for (Eigen::Index k = 0; k < energies.interCoresBond.mean().rows(); ++k)
            sum += energies.interCoresBond.mean()[k];

        for (Eigen::Index k = 0; k < energies.interCoresRest.mean().rows(); ++k)
            sum += energies.interCoresRest.mean()[k];

        if(energies.interCoresCore.mean().rows() > 0) {
            for (Eigen::Index k = 0; k < energies.interCoresCore.mean().rows() - 1; ++k) {
                for (Eigen::Index l = k + 1; l < energies.interCoresCore.mean().rows(); ++l) {
                    sum += energies.interCoresCore.mean()(k, l);
                }
            }
        }
        return sum;
    }

    void checkTotalEnergy(Group &maxima, const std::vector<Sample> &samples, std::vector<size_t> selectedNuclei, long selectedElectronsCount) const {

        auto nuclei = maxima.representative()->nuclei();
        auto repMax =  maxima.representative()->maximum();
        auto ne = repMax.numberOfEntities();
        auto na = nuclei.numberOfEntities();

        Eigen::VectorXd Te = Eigen::VectorXd::Zero(ne);
        Eigen::MatrixXd Vee = Eigen::MatrixXd::Zero(ne, ne);
        Eigen::MatrixXd Ven = Eigen::MatrixXd::Zero(ne, na);
        Eigen::MatrixXd Vnn = Eigen::MatrixXd::Zero(na, na);


        for (const auto& s : samples) {
            Te += s.kineticEnergies_;
            Vee += CoulombPotential::energies(s.sample_);
            Ven += CoulombPotential::energies(s.sample_, nuclei);
            Vnn += CoulombPotential::energies(nuclei);
        }
        // calculate mean
        Te /= samples.size();
        Vee /= samples.size();
        Ven /= samples.size();
        Vnn /= samples.size();

        maxima.sortAll();

        maxima.setSelectedElectronsCount(selectedElectronsCount);
        for(auto m : maxima){
            m.setSelectedElectronsCount(selectedElectronsCount);
        }

        LocalParticleEnergiesCalculator calculator(samples, nuclei, selectedNuclei, selectedElectronsCount);
        calculator.add(maxima);

        double eps = 1e-8;
        ASSERT_NEAR(Te.sum(), sumAll(calculator.localBondEnergies.Te), eps);
        ASSERT_NEAR(0.5 * Vee.sum(), sumAll(calculator.localBondEnergies.Vee), eps);
        ASSERT_NEAR(Ven.sum(), sumAll(calculator.localBondEnergies.Ven), eps);
        ASSERT_NEAR(0.5 * Vnn.sum(), sumAll(calculator.localBondEnergies.Vnn), eps);

        double ERef = Te.sum() + 0.5 * Vee.sum() + Ven.sum() + 0.5 * Vnn.sum();
        ASSERT_NEAR(ERef, sumAll(calculator.localBondEnergies.E), eps);
    }
};

TEST_F(ALocalParticleEnergiesCalculationTest, H2) {

    AtomsVector nuclei({
        {Element::H, {0, 0, 1}},
        {Element::H, {0, 0, -1}}});

    ElectronsVector max({
        {Spin::alpha, nuclei[0].position()},
        {Spin::beta, nuclei[1].position()}});

    ElectronsVector sample1({
        {Spin::alpha, {0, 0, 0.25}},
        {Spin::beta, {0, 0, -0.25}}});

    ElectronsVector sample2({
        {Spin::alpha, {0, 0, 0.5}},
        {Spin::beta, {0, 0, -0.5}}});

    Group maxima({
        Group({nuclei, 1.00, max, 0}),
        Group({nuclei, 1.01, max, 1})});

    std::vector<Sample> samples = {
            Sample(sample1, Eigen::VectorXd::Constant(sample1.numberOfEntities(), 0.4)),
            Sample(sample2, Eigen::VectorXd::Constant(sample2.numberOfEntities(), 0.8))
    };

    std::vector<size_t> selectedNuclei = {0, 1};
    checkTotalEnergy(maxima, samples, selectedNuclei, 2);
    checkTotalEnergy(maxima, samples, selectedNuclei, 1);
    checkTotalEnergy(maxima, samples, selectedNuclei, 0);

    selectedNuclei = {0};
    checkTotalEnergy(maxima, samples, selectedNuclei, 2);
    checkTotalEnergy(maxima, samples, selectedNuclei, 1);
    checkTotalEnergy(maxima, samples, selectedNuclei, 0);

    selectedNuclei = {};
    checkTotalEnergy(maxima, samples, selectedNuclei, 2);
    checkTotalEnergy(maxima, samples, selectedNuclei, 1);
    checkTotalEnergy(maxima, samples, selectedNuclei, 0);

}

TEST_F(ALocalParticleEnergiesCalculationTest, DISABLED_B2_not_passing) {
    AtomsVector nuclei({
                               {Element::B, {0, 0, 1}},
                               {Element::B, {0, 0, -1}}});

    ElectronsVector max({
                                  {Spin::alpha, {0, 0, 0.5}}, // => 0 and 1 are core electrons
                                  {Spin::alpha, {0, 0, -0.5}},
                                  {Spin::alpha, nuclei[0].position()},
                                  {Spin::alpha,  nuclei[1].position()}});

    ElectronsVector sample1({
                                      {Spin::alpha, {0, 0, 0.25}},
                                      {Spin::alpha, {0, 0, -0.25}},
                                      {Spin::alpha, {0, 0, 0.75}},
                                      {Spin::alpha, {0, 0, -0.75}}});

    ElectronsVector sample2 ({
                                      {Spin::alpha, {0, 0, 0.5}},
                                      {Spin::alpha, {0, 0, -0.5}},
                                      {Spin::alpha, {0, 0, 0.75}},
                                      {Spin::alpha, {0, 0, -0.75}}});


    Group maxima({
        Group({nuclei, 1.00, max, 0}),
        Group({nuclei, 1.01, max, 1})});

    std::vector<Sample> samples = {
            Sample(sample1, Eigen::VectorXd::Constant(sample1.numberOfEntities(), 0.4)),
            Sample(sample2, Eigen::VectorXd::Constant(sample2.numberOfEntities(), 0.8))
    };

    std::vector<size_t> selectedNuclei = {0};
    // the sum of interCoresBond is the exceeding amount
    checkTotalEnergy(maxima, samples, selectedNuclei, 4);

    selectedNuclei = {};
    checkTotalEnergy(maxima, samples, selectedNuclei, 3);
    checkTotalEnergy(maxima, samples, selectedNuclei, 4);

}

TEST_F(ALocalParticleEnergiesCalculationTest, B2_passing) {
    AtomsVector nuclei({
                               {Element::B, {0, 0, 1}},
                               {Element::B, {0, 0, -1}}});

    ElectronsVector max({
                                {Spin::alpha, {0, 0, 0.5}}, // => 0 and 1 are core electrons
                                {Spin::alpha, {0, 0, -0.5}},
                                {Spin::alpha, nuclei[0].position()},
                                {Spin::alpha,  nuclei[1].position()}});

    ElectronsVector sample1({
                                    {Spin::alpha, {0, 0, 0.25}},
                                    {Spin::alpha, {0, 0, -0.25}},
                                    {Spin::alpha, {0, 0, 0.75}},
                                    {Spin::alpha, {0, 0, -0.75}}});

    ElectronsVector sample2 ({
                                     {Spin::alpha, {0, 0, 0.5}},
                                     {Spin::alpha, {0, 0, -0.5}},
                                     {Spin::alpha, {0, 0, 0.75}},
                                     {Spin::alpha, {0, 0, -0.75}}});

    Group maxima({
                         Group({nuclei, 1.00, max, 0}),
                         Group({nuclei, 1.01, max, 1})});

    std::vector<Sample> samples = {
            Sample(sample1, Eigen::VectorXd::Constant(sample1.numberOfEntities(), 0.4)),
            Sample(sample2, Eigen::VectorXd::Constant(sample2.numberOfEntities(), 0.8))
    };

    std::vector<size_t> selectedNuclei = {0, 1};
    // selectedElectronsCount => the first n electrons that were permuted to the head of the vector by the local clusterer

    checkTotalEnergy(maxima, samples, selectedNuclei, 4);
    checkTotalEnergy(maxima, samples, selectedNuclei, 3);
    checkTotalEnergy(maxima, samples, selectedNuclei, 2);
    checkTotalEnergy(maxima, samples, selectedNuclei, 1);
    checkTotalEnergy(maxima, samples, selectedNuclei, 0);

    selectedNuclei = {0};
    checkTotalEnergy(maxima, samples, selectedNuclei, 3);
    checkTotalEnergy(maxima, samples, selectedNuclei, 2);
    checkTotalEnergy(maxima, samples, selectedNuclei, 1);
    checkTotalEnergy(maxima, samples, selectedNuclei, 0);


    selectedNuclei = {};
    checkTotalEnergy(maxima, samples, selectedNuclei, 2);
    checkTotalEnergy(maxima, samples, selectedNuclei, 1);
    checkTotalEnergy(maxima, samples, selectedNuclei, 0);
}