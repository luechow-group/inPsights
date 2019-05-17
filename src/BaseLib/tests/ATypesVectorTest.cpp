//
// Created by Michael Heuer on 25.04.18.
//

#include <gmock/gmock.h>
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

TEST_F(ATypesVectorTest, EigenVectorConstructor){
    SpinTypesVector stv(Eigen::VectorXi::Zero(5));
    std::stringstream stringstream;
    stringstream << stv;
    ASSERT_EQ(stringstream.str(), "-\n-\n-\n-\n-\n");
}

TEST_F(ATypesVectorTest, SpecializedConstructor){
    SpinTypesVector stv(2,3);
    std::stringstream stringstream;
    stringstream << stv;
    ASSERT_EQ(stringstream.str(), "a\na\nb\nb\nb\n");

    SpinTypesVector stv2(5);
    std::stringstream stringstream2;
    stringstream2 << stv2;
    ASSERT_EQ(stringstream2.str(), "a\na\na\na\na\n");
}

TEST_F(ATypesVectorTest,GetIndexedTypeFromIndex){
    auto indexedType = etv.getEnumeratedTypeByIndex(3);

    ASSERT_EQ(indexedType.number_,1);
    ASSERT_EQ(indexedType.type_, Element::He);
}

TEST_F(ATypesVectorTest,CheckIndexedType_Present){
    auto indexedType = EnumeratedType<Element>(Element::He,1);

    auto foundIndexPair = etv.findIndexOfEnumeratedType(indexedType);
    ASSERT_TRUE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,3);
}

TEST_F(ATypesVectorTest,CheckIndexedTypePresentBut_WrongIndex){
    auto wrongIndexedType = EnumeratedType<Element>(Element::He,3);

    auto foundIndexPair = etv.findIndexOfEnumeratedType(wrongIndexedType);
    ASSERT_FALSE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,0);
}

TEST_F(ATypesVectorTest,CheckIndexedType_MissingType){
    auto indexedType = EnumeratedType<Element>(Element::Ca,1);

    auto foundIndexPair = etv.findIndexOfEnumeratedType(indexedType);
    ASSERT_FALSE(foundIndexPair.first);
    ASSERT_EQ(foundIndexPair.second,0);
}

TEST_F(ATypesVectorTest, CountTypes_ElementTypes){
    auto result = etv.countTypes();

    ASSERT_EQ(result[0].type_, Element::H);
    ASSERT_EQ(result[0].number_, 1);

    ASSERT_EQ(result[1].type_, Element::He);
    ASSERT_EQ(result[1].number_, 3);

    ASSERT_EQ(result[2].type_, Element::Og);
    ASSERT_EQ(result[2].number_, 1);
}

TEST_F(ATypesVectorTest, CountTypes_SpinTypes){
    auto result = stv.countTypes();

    ASSERT_EQ(result[0].type_, Spin::alpha);
    ASSERT_EQ(result[0].number_, 3);

    ASSERT_EQ(result[1].type_, Spin::beta);
    ASSERT_EQ(result[1].number_, 3);
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

TEST_F(ATypesVectorTest, Multiplicity) {
    ASSERT_EQ(stv.multiplicity(),1);
    ASSERT_EQ(stvsmall.multiplicity(),2);
}

TEST_F(ATypesVectorTest, FlipSpins){
    
    auto stvFlipped = stv;
    stvFlipped.flipSpins();
    
    SpinTypesVector reference({
        Spin::beta,
        Spin::beta,
        Spin::alpha,
        Spin::beta,
        Spin::alpha,
        Spin::alpha});
    
    ASSERT_EQ(stvFlipped, reference);
}