//
// Created by Michael Heuer on 29.10.17.
//

#include <gtest/gtest.h>
#include "ParticleCollection.h"

using namespace testing;

class AParticleCollectionTest : public Test {
public:
    void SetUp() override {
        position0 = Vector3d(1,2,3);
        position1 = Vector3d(4,5,6);
        position2 = Vector3d(7,8,9);
    }
    Vector3d position0,position1,position2;
};

TEST_F(AParticleCollectionTest, IndexOperator){

    VectorXd positions(9);
    positions<< position0, position1, position2;
    ParticleCollection particleCollection(positions);

    ASSERT_EQ(particleCollection[0].position(),position0);
    ASSERT_EQ(particleCollection[1].position(),position1);
    ASSERT_EQ(particleCollection[2].position(),position2);
    ASSERT_EQ(particleCollection[-1].position(),position2);
    ASSERT_EQ(particleCollection[-2].position(),position1);
    ASSERT_EQ(particleCollection[-3].position(),position0);

    EXPECT_DEATH(particleCollection[3],"");
    EXPECT_DEATH(particleCollection[-4],"");
}

TEST_F(AParticleCollectionTest, insert){
    VectorXd positions(6);
    positions<< position0, position2;
    Particle particle(position1);
    ParticleCollection particleCollection(positions);

    particleCollection.insert(particle,1);

    ASSERT_EQ(particleCollection[0].position(),position0);
    ASSERT_EQ(particleCollection[1].position(),position1);
    ASSERT_EQ(particleCollection[2].position(),position2);
}


TEST_F(AParticleCollectionTest, prepend) {
    VectorXd positions(6);
    positions<< position1, position2;
    Particle particle(position0);
    ParticleCollection particleCollection(positions);

    particleCollection.prepend(particle);

    ASSERT_EQ(particleCollection[0].position(),position0);
    ASSERT_EQ(particleCollection[1].position(),position1);
    ASSERT_EQ(particleCollection[2].position(),position2);
}

TEST_F(AParticleCollectionTest, append) {
    VectorXd positions(6);
    positions<< position0, position1;
    Particle particle(position2);
    ParticleCollection particleCollection(positions);

    particleCollection.append(particle);

    ASSERT_EQ(particleCollection[0].position(),position0);
    ASSERT_EQ(particleCollection[1].position(),position1);
    ASSERT_EQ(particleCollection[2].position(),position2);
}


/*TEST_F(AParticleCollectionTest, extractSet){
    EXPECT_TRUE(false);
}*/