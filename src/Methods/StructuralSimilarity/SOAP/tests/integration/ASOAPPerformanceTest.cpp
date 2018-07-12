//
// Created by heuer on 02.07.18.
//

#include <gtest/gtest.h>
#include "StructuralSimilarity.h"
#include "TestMolecules.h"

class ASOAPPerformanceTest : public ::testing::Test {
public:
    double eps = std::numeric_limits<double>::epsilon()*1e3;
    double regularizationParameter = 0.5;

    void SetUp() override {
        ExpansionSettings::defaults();
        ExpansionSettings::mode = ExpansionSettings::Mode::alchemical;
    }
};

TEST_F(ASOAPPerformanceTest, Performance){

    printf("{\n");
    for (unsigned nmax = 2; nmax <= 10; ++nmax) {
        unsigned lmax = nmax;
        unsigned nParticles = 30; // nParticles = nTypes
        printf("{%d,\n",nmax);
        //calculate a molecular spectrum
        ExpansionSettings::Radial::nmax = nmax;
        ExpansionSettings::Angular::lmax = lmax;
        //ExpansionSettings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 0.5;

        //Distribute the particles on a unit circle
        double radius = 1.0;

        MolecularGeometry mol;
        for (unsigned i = 0; i < nParticles; ++i) {
            double angle = 2. * M_PI * double(i) / (nParticles - 1);
            mol.atoms().append({Element(i + 1), {cos(angle), sin(angle), 0}});
            //std::cout << mol << std::endl;
            ParticleKit::create(mol);

            double start = omp_get_wtime();
            MolecularSpectrum ms(mol);
            double t_spec = omp_get_wtime() - start;
            //printf("i=%d\nMS,elapsed time: %fs,",i,omp_get_wtime()-start);

            start = omp_get_wtime();
            double K = StructuralSimilarity::kernel(ms, ms, regularizationParameter);
            double t_sim = omp_get_wtime() - start;
            //printf(" SS=%f, elapsed time: %fs\n",K,omp_get_wtime()-start);
            printf("{%d,%f,%f},\n", i+1, t_spec, t_sim);
            fflush(stdout);
        }
        printf("},\n");
    }
    printf("}\n");
    EXPECT_TRUE(false);
}