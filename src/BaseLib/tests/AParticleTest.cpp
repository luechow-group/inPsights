//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "Particle.h"

using namespace testing;
using namespace Eigen;

class AParticleTest : public Test {
public:
    void SetUp() override {
    }
};


TEST_F(AParticleTest, Constructor){
    Particle particle1 = Particle(Vector3d(1,2,3));
    Particle particle2(Vector3d(1,2,3));
    ASSERT_TRUE(particle1.position() == particle2.position());
}

TEST_F(AParticleTest, CopyConstructor){
    Particle particle = Particle(Vector3d(1,2,3));
    Particle copyParticle(particle);
    ASSERT_TRUE(particle.position() == copyParticle.position());
}
