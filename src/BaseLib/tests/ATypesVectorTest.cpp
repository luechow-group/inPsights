//
// Created by Michael Heuer on 25.04.18.
//

#include <gtest/gtest.h>
#include <TypesVector.h>
#include <SpinType.h>
#include <ElementType.h>
#include <sstream>

using namespace testing;

class ATypesVectorTest : public Test {
public:

    SpinTypesVector stv;
    ElementTypesVector etv;

    void SetUp() override {

        stv = SpinTypesVector({Spins::SpinType::alpha,
                               Spins::SpinType::alpha,
                               Spins::SpinType::beta});
        etv = ElementTypesVector({Elements::ElementType::H,
                                  Elements::ElementType::He,
                                  Elements::ElementType::Og,
                                  Elements::ElementType::He});
    }
};

TEST_F(ATypesVectorTest, SpinTypesStream){
    std::stringstream stringstream;
    stringstream << stv;

    ASSERT_EQ(stringstream.str(), "a\na\nb\n");
}

TEST_F(ATypesVectorTest, ElementTypesStream){
    std::stringstream stringstream;
    stringstream << etv;

    ASSERT_EQ(stringstream.str(), "H\nHe\nOg\nHe\n");
}

TEST_F(ATypesVectorTest, SpecializedConstructor){
    SpinTypesVector stv(2,3);
    std::stringstream stringstream;
    stringstream << stv;
    ASSERT_EQ(stringstream.str(), "a\na\nb\nb\nb\n");
}

TEST_F(ATypesVectorTest,GetIndexedTypeFromIndex){
    auto indexedType = etv.getIndexedTypeByIndex(3);

    ASSERT_EQ(indexedType.index_,1);
    ASSERT_EQ(indexedType.type_, Elements::ElementType::He);
}

TEST_F(ATypesVectorTest,CheckIndexedType_Present){
    auto indexedType = IndexedType<Elements::ElementType>(Elements::ElementType::He,1);

    auto foundIndexPair = etv.findIndexOfIndexedType(indexedType);
    ASSERT_TRUE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,3);
}

TEST_F(ATypesVectorTest,CheckIndexedTypePresentBut_WrongIndex){
    auto wrongIndexedType = IndexedType<Elements::ElementType>(Elements::ElementType::He,2);

    auto foundIndexPair = etv.findIndexOfIndexedType(wrongIndexedType);
    ASSERT_FALSE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,0);
}

TEST_F(ATypesVectorTest,CheckIndexedType_MissingType){
    auto indexedType = IndexedType<Elements::ElementType>(Elements::ElementType::Ca,1);

    auto foundIndexPair = etv.findIndexOfIndexedType(indexedType);
    ASSERT_FALSE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,0);
}

