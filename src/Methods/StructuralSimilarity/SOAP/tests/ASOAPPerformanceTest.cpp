//
// Created by heuer on 02.07.18.
//

#include <gtest/gtest.h>
#include "StructuralSimilarity.h"
#include "TestMolecules.h"
#include <omp.h>

class ASOAPPerformanceTest : public ::testing::Test {
public:
    double eps = std::numeric_limits<double>::epsilon()*1e3;
    double regularizationParameter = 0.5;

    void SetUp() override {
        SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
    }
};

TEST_F(ASOAPPerformanceTest, DISABLED_SingleGlobalSimTiming) {

    unsigned nmax = 8;
    unsigned lmax = nmax;
    unsigned nParticles = 20; // nParticles = nTypes

    Radial::settings.nmax = nmax;
    Angular::settings.lmax = lmax;
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
    //ExpansionSettings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 0.5;


    //Distribute the particles on a unit circle
    MolecularGeometry mol;
    for (unsigned i = 0; i < nParticles; ++i) {
        double angle = 2. * M_PI * double(i) / (nParticles - 1);
        mol.atoms().append({Element(i), {cos(angle), sin(angle), 0}});
    }
    ParticleKit::create(mol);

    double start;

    start = omp_get_wtime();
    MolecularSpectrum ms(mol);
    double t1 = omp_get_wtime() - start;
    start = omp_get_wtime();
    double K = StructuralSimilarity::kernel(ms, ms, regularizationParameter);
    double t2 = omp_get_wtime() - start;

    printf("{%d,%d,%f,%f,%f},\n",nmax,nParticles,K,t1,t2);
    fflush(stdout);

    EXPECT_TRUE(false);
}


TEST_F(ASOAPPerformanceTest, DISABLED_GlobalSimPerformance){

    printf("{\n");
    for (unsigned nmax = 2; nmax <= 10; nmax=nmax+2) {
        unsigned lmax = nmax;
        unsigned nParticles = 30; // nParticles = nTypes
        printf("{%d,\n",nmax);
        //calculate a molecular spectrum
        Radial::settings.nmax = nmax;
        Angular::settings.lmax = lmax;
        //ExpansionSettings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 0.5;

        //Distribute the particles on a unit circle
        MolecularGeometry mol;
        unsigned nSkip = 5;
        for (unsigned i = 0; i < nParticles; i=i+nSkip) {
            for (unsigned j = 1; j <= nSkip; ++j) {
                double angle = 2. * M_PI * double(i+j) / (nParticles - 1);
                mol.atoms().append({Element(i + j), {cos(angle), sin(angle), 0}});
            }
            ParticleKit::create(mol);

            double start;
            MolecularSpectrum ms;

            SOAPExpansion::settings.mode = SOAPExpansion::Mode::typeAgnostic;
            start = omp_get_wtime();
            ms = MolecularSpectrum(mol);
            double t1a = omp_get_wtime() - start;
            start = omp_get_wtime();
            double generic = StructuralSimilarity::kernel(ms, ms, regularizationParameter);
            double t1b = omp_get_wtime() - start;

            SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
            start = omp_get_wtime();
            ms = MolecularSpectrum(mol);
            double t2a = omp_get_wtime() - start;
            start = omp_get_wtime();
            double chemical = StructuralSimilarity::kernel(ms, ms, regularizationParameter);
            double t2b = omp_get_wtime() - start;

            SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
            start = omp_get_wtime();
            ms = MolecularSpectrum(mol);
            double t3a = omp_get_wtime() - start;
            start = omp_get_wtime();
            double alchemical = StructuralSimilarity::kernel(ms, ms, regularizationParameter);
            double t3b = omp_get_wtime() - start;

            printf("{%d,%f,%f,%f,%f,%f,%f},\n", i+nSkip,t1a,t2a,t3a,t1b,t2b,t3b);
            fflush(stdout);
        }
        printf("},\n");
    }
    printf("}\n");
    EXPECT_TRUE(false);
}

TEST_F(ASOAPPerformanceTest, DISABLED_LocalSimPerformance){

    printf("{\n");
    for (unsigned nmax = 2; nmax <= 10; nmax=nmax+2) {
        unsigned lmax = nmax;
        unsigned nParticles = 30; // nParticles = nTypes
        printf("{%d,\n",nmax);
        //calculate a molecular spectrum
        Radial::settings.nmax = nmax;
        Angular::settings.lmax = lmax;
        //ExpansionSettings::Alchemical::pairSimilarities[{int(Spin::alpha),int(Spin::beta)}] = 0.5;

        //Distribute the particles on a unit circle
        double radius = 1.0;

        MolecularGeometry mol;
        unsigned nSkip = 5;
        for (unsigned i = 0; i < nParticles; i=i+nSkip) {
            for (unsigned j = 1; j <= nSkip; ++j) {
                double angle = 2. * M_PI * double(i+j) / (nParticles - 1);
                mol.atoms().append({Element(i + j), {cos(angle), sin(angle), 0}});
            }
            ParticleKit::create(mol);


            double start;
            Environment e1{{},{}};

            MolecularSpectrum ms(mol);

            SOAPExpansion::settings.mode = SOAPExpansion::Mode::typeAgnostic;
            e1 = Environment(mol,mol.atoms()[0].position());
            start = omp_get_wtime();
            double generic = LocalSimilarity::kernel(e1,e1);
            double t1 = omp_get_wtime() - start;

            SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;
            e1 = Environment(mol,mol.atoms()[0].position());
            start = omp_get_wtime();
            double chemical = LocalSimilarity::kernel(e1,e1);
            double t2 = omp_get_wtime() - start;

            SOAPExpansion::settings.mode = SOAPExpansion::Mode::alchemical;
            e1 = Environment(mol,mol.atoms()[0].position());
            start = omp_get_wtime();
            double alchemical = LocalSimilarity::kernel(e1,e1);
            double t3 = omp_get_wtime() - start;

            //printf(" SS=%f, elapsed time: %fs\n",K,omp_get_wtime()-start);
            printf("{%d,%f,%f,%f},\n", i+nSkip,t1,t2,t3);
            fflush(stdout);
        }
        printf("},\n");
    }
    printf("}\n");
    EXPECT_TRUE(false);
}