//
// Created by Michael Heuer on 25.04.18.
//

#include <gtest/gtest.h>
#include <TypesVector.h>
#include <SpinType.h>
#include <ElementType.h>
#include <sstream>

using namespace testing;
using namespace Spins;

class ATypesVectorTest : public Test {
public:
    void SetUp() override {
    }
};

TEST_F(ATypesVectorTest, SpinTypesStream){
    TypesVector<Spins::SpinType> tv({Spins::SpinType::alpha,
                                     Spins::SpinType::alpha,
                                     Spins::SpinType::beta});
    std::stringstream stringstream;
    stringstream << tv;

    ASSERT_TRUE(stringstream.str() == "a\na\nb\n");
}

TEST_F(ATypesVectorTest, ElementTypesStream){
    TypesVector<Elements::ElementType> tv({Elements::ElementType::H,
                                           Elements::ElementType::He,
                                           Elements::ElementType::Og});
    std::stringstream stringstream;
    stringstream << tv;
    std::cout << tv;

    ASSERT_TRUE(stringstream.str() == "H\nHe\nOg\n");
}

TEST_F(ATypesVectorTest, SpecializedConstructor){
    SpinTypesVector stv(2,3);
    std::stringstream stringstream;
    stringstream << stv;
    std::cout << stv;

}