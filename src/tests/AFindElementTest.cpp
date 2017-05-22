//
// Created by Moria on 22.05.2017.
//

#include "pse.h"
#include "gtest/gtest.h"

TEST(FindElementTest, ExistingElements)
{
    std::string testElements[11]= {"C", "H", "Xe", "Pb", "Kr", "Ni", "Cl", "At", "Ra", "Lu", "Hg"};
    int testElementOZ[11]={6,1,54,82,36,28,17,85,88,71,80};
    for(int i=0;i<11;i++){
        EXPECT_EQ(testElementOZ[i],Pse::findElement(testElements[i]));
    }
}

TEST(FindElementTest, NotExistingElements)
{
    std::string testElements[3]= {"L", "Osdorf", "Ze"};
    for(int i=0;i<3;i++){
        EXPECT_EQ(0,Pse::findElement(testElements[i]));
    }
}



