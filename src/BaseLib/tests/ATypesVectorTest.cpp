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

/*TEST_F(ATypesVectorTest, CopyConstructor){//TODO
    SpinTypesVector s(stv);
    bool b = (s==stv);
    std::cout << b << std::endl;
}*/

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

TEST_F(ATypesVectorTest, Slice){
    auto s = stvsmall;
    auto types = s.typesAsEigenVector();

    ASSERT_EQ(s.entity(0).typesRef(), types.segment(0,1));
    ASSERT_EQ(s.entity(1).typesRef(), types.segment(1,1));
    ASSERT_EQ(s.entity(2).typesRef(), types.segment(2,1));

    ASSERT_EQ(s.slice({0,1}).typesRef(),s.entity(0).typesRef());
    ASSERT_EQ(s.slice({1,1}).typesRef(),s.entity(1).typesRef());
    ASSERT_EQ(s.slice({2,1}).typesRef(),s.entity(2).typesRef());

    ASSERT_EQ(s.slice({0,2}).typesRef(), types.segment(0,2));
    ASSERT_EQ(s.slice({1,2}).typesRef(), types.segment(1,2));
    ASSERT_EQ(s.slice({0,3}).typesRef(), types.segment(0,3));

    ASSERT_EQ(s.slice({0,3}).typesRef(),s.typesRef());

    EXPECT_DEATH(s.entity(-1).typesRef(),"");
    EXPECT_DEATH(s.slice({0,4}).typesRef(),"");
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

TEST_F(ATypesVectorTest, PermuteSlice){
    auto s = stvsmall;

    VectorXi p(2);
    p << 1,0;

    s.slice({1,2}).permute(PermutationMatrix<Dynamic>(p));
    ASSERT_EQ(s[0],Spin::alpha);
    ASSERT_EQ(s[1],Spin::beta);
    ASSERT_EQ(s[2],Spin::alpha);

}

TEST_F(ATypesVectorTest, EqualityOperator) {
    SpinTypesVector s1({Spin::alpha, Spin::alpha, Spin::beta});
    SpinTypesVector s2({Spin::alpha, Spin::alpha, Spin::beta});
    SpinTypesVector s3({Spin::beta, Spin::alpha, Spin::beta});

    ASSERT_TRUE(s1 == s2);
    ASSERT_TRUE(s1.slice({0, 2}) == s2.slice({0, 2}));

    ASSERT_FALSE(s1.slice({0, 2}) == s2.slice({1, 2}));
    ASSERT_FALSE(s1 == s3);

    //equality means that not only the slice but also the underlying data is identical
    ASSERT_FALSE(s1.slice({1, 2}) == s3.slice({1, 2}));
}