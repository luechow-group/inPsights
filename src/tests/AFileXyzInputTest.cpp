//
// Created by Moria on 22.05.2017.
//

#include <HungarianElectronAssigner.h>
#include "gmock/gmock.h"
#include "Molecule.h"
#include "FileXyzInput.h"
//#include "Eigen/dense"

TEST(TestFileXyzInput,MoleculeCores)
{
    Molecule newMolecule;
    FileXyzInput input("../../../testinput/EPA.ref","../../../testinput/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    std::vector<Eigen::Vector3d> expectedPositions;
    expectedPositions.emplace_back();
    (*expectedPositions.rbegin())<<0,0,0;
    expectedPositions.emplace_back();
    (*expectedPositions.rbegin())<<1.1892424,-1.1892235,-1.1892235;
    expectedPositions.emplace_back();
    (*expectedPositions.rbegin())<<-1.1892235,1.1892046,-1.1892235;
    expectedPositions.emplace_back();
    (*expectedPositions.rbegin())<<-1.1892235,-1.1892235,1.1892235;
    expectedPositions.emplace_back();
    (*expectedPositions.rbegin())<<1.1892235,1.1892235,1.1892235;
    std::vector<std::string> expectedElements={"C","H","H","H","H"};
    std::vector<int> expectedOZ={6,1,1,1,1};
    EXPECT_EQ(newMolecule.getCores().size(),5);
    for(int i=0;i<5;i++){
        EXPECT_TRUE(expectedPositions[i].isApprox(newMolecule.getCores()[i].getPosition()));
        EXPECT_EQ(expectedOZ[i],newMolecule.getCores()[i].getCharge());
        EXPECT_EQ(expectedElements[i],newMolecule.getCores()[i].getElementType());
    }
}

TEST(TestFileXyzInput, FileNameNotExistent)
{
    Molecule newMolecule;
    EXPECT_EXIT(FileXyzInput input("ulfbaggabagga","ulfbaggabgaa"),::testing::ExitedWithCode(1), "does not exist");
}

TEST(TestFileXyzInput, ReadAssignment)
{
    Molecule newMolecule;
    FileXyzInput input("../../../testinput/EPA.ref","../../../testinput/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner hea;
    input.readElectronCoreAssignations(newMolecule.getCores(),hea);
    EXPECT_EQ(input.getAssignations().size(),15);
    EXPECT_THAT(input.getAssignations()[0][0].second,::testing::UnorderedElementsAre(0,1,2,3,4,5));
    EXPECT_THAT(input.getAssignations()[0][1].second,::testing::UnorderedElementsAre(6));
    EXPECT_THAT(input.getAssignations()[0][2].second,::testing::UnorderedElementsAre(7));
    EXPECT_THAT(input.getAssignations()[0][3].second,::testing::UnorderedElementsAre(8));
    EXPECT_THAT(input.getAssignations()[0][4].second,::testing::UnorderedElementsAre(9));
}
