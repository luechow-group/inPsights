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
            {{Elements::ElementType::H,{0,0, 0.37}},
             {Elements::ElementType::H,{0,0,-0.37}}}),
         ElectronsVector(
                 {{Spins::SpinType::alpha,{0,0, 0.37}},
                  {Spins::SpinType::alpha,{0,0,-0.37}}})
    };
    B = {AtomsVector(
            {{Elements::ElementType::H,{0,0, 0.37}},
             {Elements::ElementType::He,{0,0,-0.37}}}),
         ElectronsVector(
                 {{Spins::SpinType::alpha,{0,0, 0.37}},
                  {Spins::SpinType::alpha,{0,0,-0.37}}})
    };

    ExpansionSettings::defaults();
    ExpansionSettings::mode = ExpansionMode::TypeSpecific;
    ParticleKit::create({{Elements::ElementType::H,2},{Elements::ElementType::He,2}},{2,2});

    double result = StructuralSimilarity::stucturalSimilarity(A,B,1);

    end = std::chrono::system_clock::now();
    std::cout << "FINAL RESULT:" << result << std::endl;

    long elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
            (end-start).count();
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds << "ms\n";
}