//
// Created by Michael Heuer on 06.09.18.
//

#include <gmock/gmock.h>
#include <HungarianHelper.h>
#include <TestMolecules.h>
#include <limits>

TEST(AHungarianHelperTest, CombinePermutations){
    Eigen::PermutationMatrix<Eigen::Dynamic>p1,p2;
    p1.setIdentity(4);
    p2.setIdentity(4);

    Eigen::VectorXi expected(8);
    expected << 0,1,2,3,4,5,6,7;

    ASSERT_TRUE(p1.indices().base().isApprox(expected.segment(0,4)));
    ASSERT_TRUE(p2.indices().base().isApprox(expected.segment(0,4)));

    auto p = HungarianHelper::combinePermutations(p1,p2);
    ASSERT_TRUE(p.indices().base().isApprox(expected));
}

TEST(AHungarianHelperTest, SpinSpecificHungarian){
    auto eNormal = TestMolecules::eightElectrons::square.electrons();

    Eigen::VectorXi alphaBetaPositionFlip(8);
    alphaBetaPositionFlip << 4,5,6,7,0,1,2,3;
    Eigen::PermutationMatrix<Eigen::Dynamic> p(alphaBetaPositionFlip);

    auto eFlipped = eNormal;
    eFlipped.positionsVector().permute(p);

    //auto bestMatch = HungarianHelper::spinSpecificHungarian(eNormal, eFlipped, false);
    auto bestMatchFlipped = HungarianHelper::spinSpecificBestMatch(eNormal, eFlipped, true);
    auto bestMatchFlippedInverse = Eigen::PermutationMatrix<Eigen::Dynamic>(bestMatchFlipped.inverse());

    ASSERT_TRUE(bestMatchFlippedInverse.indices().base().isApprox(p.indices().base()));
}

TEST(AHungarianHelperTest, BestMatchNorm) {
    ElectronsVector v1({
        {Spin::alpha, {0,1,2}},
        {Spin::alpha, {0,0,0}}});

    ElectronsVector v2({
            {Spin::alpha, {0,0,0}},
            {Spin::alpha, {0,4,6}}});


    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 1,0;

    auto [norm2, perm2] = Metrics::bestMatch<2, 2>(v1, v2);
    ASSERT_EQ(norm2, 5.0);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto [normInf, permInf] = Metrics::bestMatch<Eigen::Infinity, 2>(v1, v2);
    ASSERT_EQ(normInf,5.0);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto [normInfInf, permInfInf] = Metrics::bestMatch<Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 4.0);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));

}

TEST(AHungarianHelperTest, SpinSpecificBestMatchNormSameSpin) {
    ElectronsVector v1({
        {Spin::alpha, {0,1,2}},
        {Spin::alpha, {0,0,0}}});

    ElectronsVector v2({
        {Spin::alpha, {0,0,0}},
        {Spin::alpha, {0,4,6}}});


    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 1,0;

    auto [norm2, perm2] = Metrics::spinSpecificBestMatch<2, 2>(v1, v2);
    ASSERT_EQ(norm2, 5.0);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto [normInf, permInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, 2>(v1, v2);
    ASSERT_EQ(normInf,5.0);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto [normInfInf, permInfInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 4.0);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));
}

TEST(AHungarianHelperTest, SpinSpecificBestMatchNormDifferentSpin) {
    ElectronsVector v1({
        {Spin::alpha, {0,1,2}},
        {Spin::beta, {0,0,0}}});

    ElectronsVector v2({
        {Spin::alpha, {0,0,0}},
        {Spin::beta, {0,4,6}}});

    Eigen::VectorXi expectedPerm(2);
    expectedPerm << 0,1;

    auto eps = std::numeric_limits<double>::epsilon()*10;

    auto [norm2, perm2] = Metrics::spinSpecificBestMatch<2, 2>(v1, v2);
    ASSERT_NEAR(norm2, std::sqrt(1+2*2+4*4+6*6), eps);
    ASSERT_TRUE(perm2.indices().isApprox(expectedPerm));

    auto [normInf, permInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, 2>(v1, v2);
    ASSERT_NEAR(normInf, std::sqrt(4*4+6*6), eps);
    ASSERT_TRUE(permInf.indices().isApprox(expectedPerm));

    auto [normInfInf, permInfInf] = Metrics::spinSpecificBestMatch<Eigen::Infinity, Eigen::Infinity>(v1, v2);
    ASSERT_EQ(normInfInf, 6);
    ASSERT_TRUE(permInfInf.indices().isApprox(expectedPerm));
}

TEST(AHungarianHelperTest, RealMaxima){
    ElectronsVector v1({
    {Spin::alpha,{-1.924799,-0.000888,-2.199093}},
    {Spin::alpha,{-0.425365,-0.687079, 1.739155}},
    {Spin::alpha,{ 0.807722, 0.024079, 1.739155}},
    {Spin::alpha,{ 0.000000, 0.000000,-1.446226}},
    {Spin::alpha,{ 0.397920,-0.688486,-1.772321}},
    {Spin::alpha,{-0.961644, 1.667362, 2.199093}},
    {Spin::alpha,{-0.019931, 0.034495,-0.660905}},
    {Spin::alpha,{ 0.000000, 0.000000, 1.446226}},
    {Spin::alpha,{ 0.961644, 1.667362,-2.199093}},
    {Spin::beta ,{ 0.000000, 0.000000,-1.446226}},
    {Spin::beta ,{ 0.000000, 0.000000, 1.446226}},
    {Spin::beta ,{-0.963174,-1.666493, 2.199093}},
    {Spin::beta ,{-0.397172, 0.688642, 1.772829}},
    {Spin::beta ,{ 0.963174,-1.666493,-2.199093}},
    {Spin::beta ,{-0.807693,-0.025002,-1.739139}},
    {Spin::beta ,{ 1.924799,-0.000888, 2.199093}},
    {Spin::beta ,{ 0.424858, 0.687373,-1.739140}},
    {Spin::beta ,{ 0.020247,-0.035100, 0.660925}}});

    ElectronsVector v2({
    {Spin::alpha,{ 0.019994, 0.034636, 0.660859}},
    {Spin::alpha,{ 0.000000, 0.000000, 1.446226}},
    {Spin::alpha,{-0.397902,-0.688437, 1.772419}},
    {Spin::alpha,{ 0.000000, 0.000000,-1.446226}},
    {Spin::alpha,{ 0.961644, 1.667362,-2.199093}},
    {Spin::alpha,{-0.961644, 1.667362, 2.199093}},
    {Spin::alpha,{ 0.425527,-0.687038,-1.739081}},
    {Spin::alpha,{-0.807764, 0.024240,-1.739086}},
    {Spin::alpha,{ 1.924799,-0.000888, 2.199093}},
    {Spin::beta ,{-0.963174,-1.666493, 2.199093}},
    {Spin::beta ,{ 0.000000, 0.000000, 1.446226}},
    {Spin::beta ,{-0.424878, 0.687387, 1.739079}},
    {Spin::beta ,{ 0.397266, 0.688792,-1.772418}},
    {Spin::beta ,{ 0.000000, 0.000000,-1.446226}},
    {Spin::beta ,{ 0.963174,-1.666493,-2.199093}},
    {Spin::beta ,{ 0.807718,-0.025010, 1.739076}},
    {Spin::beta ,{-1.924799,-0.000888,-2.199093}},
    {Spin::beta ,{-0.019967,-0.034674,-0.660818}}});

    auto [norm, perm] = Metrics::bestMatch<Eigen::Infinity, 2>(v1, v2);

    Eigen::VectorXi expectedPerm(v1.numberOfEntities());
    expectedPerm << 16,2,15,3,6,5,17,1,4,13,10,9,11,14,7,8,12,0;

    ASSERT_LT(norm, 0.1);
    ASSERT_TRUE(perm.indices().isApprox(expectedPerm));
}


TEST(AHungarianHelperTest, BestMatchSimilarity_threeElectrons){
    MolecularGeometry A = {
            AtomsVector(),
            ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                             {Spin::alpha, {0, 0, 0.5}},
                             {Spin::beta,  {0, 0, 0.3}}
            })};
    MolecularGeometry B = {
            AtomsVector(),
            ElectronsVector({{Spin::beta,  {0, 0, 0.3}},
                             {Spin::alpha, {0, 0, 0.5}},
                             {Spin::alpha, {0, 0,-0.5}}
            })};
    MolecularGeometry C = {
            AtomsVector(),
            ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                             {Spin::alpha, {0, 0, 0.5}},
                             {Spin::beta,  {0, 0,-0.3}}
            })};
    MolecularGeometry D = {
            AtomsVector(),
            ElectronsVector({{Spin::alpha, {0, 0,-0.5}},
                             {Spin::beta,  {0, 0,-0.3}},
                             {Spin::alpha, {0, 0, 0.5}}
            })};

    ParticleKit::create(A);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);
    auto specC = MolecularSpectrum(C);
    auto specD = MolecularSpectrum(D);

    auto idPerm = Eigen::PermutationMatrix<Eigen::Dynamic>(A.electrons().numberOfEntities());
    idPerm.setIdentity();


    auto [normAtoA, permAtoA] = Metrics::bestMatchSimilarity(specA, specA);
    ASSERT_EQ(normAtoA,1);
    ASSERT_TRUE(permAtoA.indices().isApprox(idPerm.indices()));

    auto [normBtoA, permBtoA] = Metrics::bestMatchSimilarity(specB, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoA(Eigen::Vector3i(2,1,0));
    ASSERT_EQ(normBtoA,1);
    ASSERT_TRUE(permBtoA.indices().isApprox(refPermBtoA.indices()));

    auto [normCtoA, permCtoA] = Metrics::bestMatchSimilarity(specC, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoA(Eigen::Vector3i(1,0,2));
    ASSERT_EQ(normCtoA,1);
    ASSERT_TRUE(permCtoA.indices().isApprox(refPermCtoA.indices()));

    auto [normDtoA, permDtoA] = Metrics::bestMatchSimilarity(specD, specA);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoA(Eigen::Vector3i(1,2,0));
    ASSERT_EQ(normDtoA,1);
    ASSERT_TRUE(permDtoA.indices().isApprox(refPermDtoA.indices()));


    auto [normAtoB, permAtoB] = Metrics::bestMatchSimilarity(specA, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoB(Eigen::Vector3i(2,1,0));
    ASSERT_EQ(normAtoB,1);
    ASSERT_TRUE(permAtoB.indices().isApprox(refPermAtoB.indices()));

    auto [normBtoB, permBtoB] = Metrics::bestMatchSimilarity(specB, specB);
    ASSERT_EQ(normBtoB,1);
    ASSERT_TRUE(permBtoB.indices().isApprox(idPerm.indices()));

    auto [normCtoB, permCtoB] = Metrics::bestMatchSimilarity(specC, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoB(Eigen::Vector3i(1,2,0));
    ASSERT_EQ(normCtoB,1);
    ASSERT_TRUE(permCtoB.indices().isApprox(refPermCtoB.indices()));

    auto [normDtoB, permDtoB] = Metrics::bestMatchSimilarity(specD, specB);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoB(Eigen::Vector3i(1,0,2));
    ASSERT_EQ(normDtoB,1);
    ASSERT_TRUE(permDtoB.indices().isApprox(refPermDtoB.indices()));


    auto [normAtoC, permAtoC] = Metrics::bestMatchSimilarity(specA, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoC(Eigen::Vector3i(1,0,2));
    ASSERT_EQ(normAtoC,1);
    ASSERT_TRUE(permAtoC.indices().isApprox(refPermAtoC.indices()));

    auto [normBtoC, permBtoC] = Metrics::bestMatchSimilarity(specB, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoC(Eigen::Vector3i(2,0,1));
    ASSERT_EQ(normBtoC,1);
    ASSERT_TRUE(permBtoC.indices().isApprox(refPermBtoC.indices()));

    auto [normCtoC, permCtoC] = Metrics::bestMatchSimilarity(specC, specC);
    ASSERT_EQ(normCtoC,1);
    ASSERT_TRUE(permCtoC.indices().isApprox(idPerm.indices()));

    auto [normDtoC, permDtoC] = Metrics::bestMatchSimilarity(specD, specC);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermDtoC(Eigen::Vector3i(0,2,1));
    ASSERT_EQ(normDtoC,1);
    ASSERT_TRUE(permDtoC.indices().isApprox(refPermDtoC.indices()));


    auto [normAtoD, permAtoD] = Metrics::bestMatchSimilarity(specA, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAtoD(Eigen::Vector3i(2,0,1));
    ASSERT_EQ(normAtoD,1);
    ASSERT_TRUE(permAtoD.indices().isApprox(refPermAtoD.indices()));

    auto [normBtoD, permBtoD] = Metrics::bestMatchSimilarity(specB, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermBtoD(Eigen::Vector3i(1,0,2));
    ASSERT_EQ(normBtoD,1);
    ASSERT_TRUE(permBtoD.indices().isApprox(refPermBtoD.indices()));

    auto [normCtoD, permCtoD] = Metrics::bestMatchSimilarity(specC, specD);
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermCtoD(Eigen::Vector3i(0,2,1));
    ASSERT_EQ(normCtoD,1);
    ASSERT_TRUE(permCtoD.indices().isApprox(refPermCtoD.indices()));

    auto [normDtoD, permDtoD] = Metrics::bestMatchSimilarity(specD, specD);
    ASSERT_EQ(normDtoD,1);
    ASSERT_TRUE(permDtoD.indices().isApprox(idPerm.indices()));
}

TEST(AHungarianHelperTest, BestMatchSimilarity_BH3){

    auto A = TestMolecules::BH3::ionicMirrored;
    auto B = TestMolecules::BH3::ionic;

    ParticleKit::create(A);
    SOAPExpansion::settings.mode = SOAPExpansion::Mode::chemical;

    auto specA = MolecularSpectrum(A);
    auto specB = MolecularSpectrum(B);

    auto idPerm = Eigen::PermutationMatrix<Eigen::Dynamic>(A.electrons().numberOfEntities());
    idPerm.setIdentity();

     auto [normAA, permAA] = Metrics::bestMatchSimilarity(specB, specB);
    ASSERT_EQ(normAA,1);
    ASSERT_TRUE(permAA.indices().isApprox(idPerm.indices()));

    auto [normBB, permBB] = Metrics::bestMatchSimilarity(specB, specB);
    ASSERT_EQ(normBB,1);
    ASSERT_TRUE(permBB.indices().isApprox(idPerm.indices()));

    auto [normAB, permAB] = Metrics::bestMatchSimilarity(specA, specB);
    std::cout << normAB << std::endl;
    std::cout << permAB.indices().transpose() << std::endl;
    Eigen::VectorXi indices(A.electrons().numberOfEntities());
    indices << 0, 1, 4, 5, 2, 3, 6, 7;
    Eigen::PermutationMatrix<Eigen::Dynamic> refPermAB(indices);
    ASSERT_EQ(normAB,1);
    ASSERT_TRUE(permAB.indices().isApprox(refPermAB.indices()));
}