//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
#include "ElectronAssigner.h"

/*
 * ATTENTION: Do not use this ElectronAssigner if you intend to get scientific restults.
 *
 * This is an Implementation of an ElectronAssigner, which is intended to use only for functionality testing.
 * The Electrons are assigned to the cores by "first comes, first served"-principle.
 *
 */
class TestElectronAssigner : public ElectronAssigner{
public:
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) override;
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons) override;
};


#endif //LOCALSPINMULTIPLICITY_TESTELECTRONASSIGNER_H
