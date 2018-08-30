//
// Created by Michael Heuer on 25.04.18.
//

#include <gtest/gtest.h>
#include <TypesVector.h>
#include <SpinType.h>
#include <ElementType.h>
#include <sstream>

using namespace testing;
using namespace Eigen;

class ATypesVectorTest : public Test {
public:

    SpinTypesVector stv, stvsmall;
    ElementTypesVector etv;

    void SetUp() override {

        stv = SpinTypesVector({Spin::alpha,
                               Spin::alpha,
                               Spin::beta,
                               Spin::alpha,
                               Spin::beta,
                               Spin::beta});

        stvsmall = SpinTypesVector({Spin::alpha,
                                    Spin::alpha,
                                    Spin::beta});

        etv = ElementTypesVector({Element::H,
                                  Element::He,
                                  Element::Og,
                                  Element::He,
                                  Element::He});
    }
};

TEST_F(ATypesVectorTest, SpinTypesStream){
    std::stringstream stringstream;
    stringstream << stv;

    ASSERT_EQ(stringstream.str(), "a\na\nb\na\nb\nb\n");
}

TEST_F(ATypesVectorTest, ElementTypesStream){
    std::stringstream stringstream;
    stringstream << etv;

    ASSERT_EQ(stringstream.str(), "H\nHe\nOg\nHe\nHe\n");
}

TEST_F(ATypesVectorTest, SpecializedConstructor){
    SpinTypesVector stv(2,3);
    std::stringstream stringstream;
    stringstream << stv;
    ASSERT_EQ(stringstream.str(), "a\na\nb\nb\nb\n");
}

TEST_F(ATypesVectorTest,GetIndexedTypeFromIndex){
    auto indexedType = etv.getNumberedTypeByIndex(3);

    ASSERT_EQ(indexedType.number_,1);
    ASSERT_EQ(indexedType.type_, Element::He);
}

TEST_F(ATypesVectorTest,CheckIndexedType_Present){
    auto indexedType = NumberedType<Element>(Element::He,1);

    auto foundIndexPair = etv.findIndexOfNumberedType(indexedType);
    ASSERT_TRUE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,3);
}

TEST_F(ATypesVectorTest,CheckIndexedTypePresentBut_WrongIndex){
    auto wrongIndexedType = NumberedType<Element>(Element::He,3);

    auto foundIndexPair = etv.findIndexOfNumberedType(wrongIndexedType);
    ASSERT_FALSE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,0);
}

TEST_F(ATypesVectorTest,CheckIndexedType_MissingType){
    auto indexedType = NumberedType<Element>(Element::Ca,1);

    auto foundIndexPair = etv.findIndexOfNumberedType(indexedType);
    ASSERT_FALSE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,0);
}

TEST_F(ATypesVectorTest, CountTypes_ElementTypes){
    auto result = etv.countTypes();

    ASSERT_EQ(result[0].first, Element::H);
    ASSERT_EQ(result[0].second, 1);

    ASSERT_EQ(result[1].first, Element::He);
    ASSERT_EQ(result[1].second, 3);

    ASSERT_EQ(result[2].first, Element::Og);
    ASSERT_EQ(result[2].second, 1);
}

TEST_F(ATypesVectorTest, CountTypes_SpinTypes){
    auto result = stv.countTypes();

    ASSERT_EQ(result[0].first, Spin::alpha);
    ASSERT_EQ(result[0].second, 3);

    ASSERT_EQ(result[1].first, Spin::beta);
    ASSERT_EQ(result[1].second, 3);
}

TEST_F(ATypesVectorTest, Permute){
    auto s = stvsmall;

    VectorXi p(3);
    p << 2,0,1;

    s.permute(PermutationMatrix<Dynamic>(p));
    ASSERT_EQ(s[0],Spin::alpha);
    ASSERT_EQ(s[1],Spin::beta);
    ASSERT_EQ(s[2],Spin::alpha);

}