//
// Created by Michael Heuer on 25.05.18.
//
#include <MolecularGeometry.h>
#include <StructuralSimilarity.h>
#include <chrono>
#include <ctime>


int main(int argc, char *argv[]) {

    std::chrono::time_point<std::chrono::system_clock> start, end;

    start = std::chrono::system_clock::now();

    MolecularGeometry A,B;
    A = {AtomsVector(
            {{Element::H,{0,0, 0.37}},
             {Element::H,{0,0,-0.37}}}),
         ElectronsVector(
                 {{Spin::alpha,{0,0, 0.37}},
                  {Spin::alpha,{0,0,-0.37}}})
    };
    B = {AtomsVector(
            {{Element::H,{0,0, 0.37}},
             {Element::He,{0,0,-0.37}}}),
         ElectronsVector(
                 {{Spin::alpha,{0,0, 0.37}},
                  {Spin::alpha,{0,0,-0.37}}})
    };

    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionSettings::Mode::TypeSpecific;
    ParticleKit::create({{Element::H,2},{Element::He,2}},{2,2});

    double result = StructuralSimilarity::stucturalSimilarity(A,B,1);

    end = std::chrono::system_clock::now();
    std::cout << "FINAL RESULT:" << result << std::endl;

    long elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
            (end-start).count();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds << "ms\n";
}