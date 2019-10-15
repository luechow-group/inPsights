/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <gmock/gmock.h>
#include <Combinatorics.h>
#include <vector>
#include <deque>

struct Thing{
    double number;
    std::string name;
};

using namespace testing;


TEST(ACombinatoricsTest, CombinationsCount) {
    ASSERT_EQ(Combinatorics::binomial(0, 0),1); //
    ASSERT_EQ(Combinatorics::binomial(0, 1),0); // {}

    ASSERT_EQ(Combinatorics::binomial(1, 0),1); // {0}
    ASSERT_EQ(Combinatorics::binomial(1, 1),1); // 0
    ASSERT_EQ(Combinatorics::binomial(1, 2),0); //

    ASSERT_EQ(Combinatorics::binomial(2, 0),1); // {01}
    ASSERT_EQ(Combinatorics::binomial(2, 1),2); // 0 1
    ASSERT_EQ(Combinatorics::binomial(2, 2),1); // 01
    ASSERT_EQ(Combinatorics::binomial(2, 3),0); //

    ASSERT_EQ(Combinatorics::binomial(3, 0),1); // {012}
    ASSERT_EQ(Combinatorics::binomial(3, 1),3); // 0 1 2
    ASSERT_EQ(Combinatorics::binomial(3, 2),3); // 01 02 12
    ASSERT_EQ(Combinatorics::binomial(3, 3),1); // 012
    ASSERT_EQ(Combinatorics::binomial(3, 4),0); //

    ASSERT_EQ(Combinatorics::binomial(4, 0),1); // {0123}
    ASSERT_EQ(Combinatorics::binomial(4, 1),4); // 0 1 2 3
    ASSERT_EQ(Combinatorics::binomial(4, 2),6); // 01 02 03 12 13 23
    ASSERT_EQ(Combinatorics::binomial(4, 3),4); // 012 013 023 123
    ASSERT_EQ(Combinatorics::binomial(4, 4),1); // 0123
    ASSERT_EQ(Combinatorics::binomial(4, 5),0); //

    ASSERT_EQ(Combinatorics::binomial(5, 0),1); // {01234}
    ASSERT_EQ(Combinatorics::binomial(5, 1),5); // 0 1 2 3 4
    ASSERT_EQ(Combinatorics::binomial(5, 2),10);// 01 02 03 04 12 13 13 23 24 34
    ASSERT_EQ(Combinatorics::binomial(5, 3),10);// 012 013 014 023 024 034 123 124 134 234
    ASSERT_EQ(Combinatorics::binomial(5, 4),5); // 0123 0124 0134 0234 1234
    ASSERT_EQ(Combinatorics::binomial(5, 5),1); // 01234
    ASSERT_EQ(Combinatorics::binomial(5, 6),0); //
}

TEST(ACombinatoricsTest, LengthZeroCombination) {
    std::vector<int> s{0,1,2,3};

    Combinatorics::Combinations<int> combinations(s,0);
    ASSERT_EQ(combinations.get().size(), 1);
    ASSERT_EQ(combinations.get()[0].size(), 0);
}

TEST(ACombinatoricsTest, IntegerCombinations) {
    std::vector<int> s{0,1,2,3};

    Combinatorics::Combinations<int> combinations(s,3);

    ASSERT_THAT(combinations.get()[0], ElementsAre(0,1,2));
    ASSERT_THAT(combinations.get()[1], ElementsAre(0,1,3));
    ASSERT_THAT(combinations.get()[2], ElementsAre(0,2,3));
    ASSERT_THAT(combinations.get()[3], ElementsAre(1,2,3));
}

TEST(ACombinatoricsTest, ThingCombinations) {
    std::vector<Thing> things =  {{1.2, "Joe"},{-2.6, "Bart"}};

    Combinatorics::Combinations<Thing> combinations(things,2);

    ASSERT_EQ(combinations.get()[0][0].number, 1.2);
    ASSERT_EQ(combinations.get()[0][1].number, -2.6);
}

TEST(ACombinatoricsTest, OutputContainer) {
    std::vector<Thing> things =  {{1.2, "Joe"},{-2.6, "Bart"}};
    Combinatorics::Combinations<Thing, std::vector<Thing>, std::deque<std::vector<Thing>>> combinations(things,1);

    ASSERT_EQ(combinations.get()[0][0].number, 1.2);
    ASSERT_EQ(combinations.get()[1][0].number, -2.6);

    combinations.get().pop_front();

    ASSERT_EQ(combinations.get()[0][0].number, -2.6);
}

TEST(ACombinatoricsTest, LengthZeroPermutation) {
    std::vector<int> s{};
    Combinatorics::Permutations<int> permutations(s);
    ASSERT_EQ(permutations.get().size(),Combinatorics::factorial(s.size()));
}


TEST(ACombinatoricsTest, LengthOnePermutation) {
    std::vector<int> s{0};
    Combinatorics::Permutations<int> permutations(s);
    ASSERT_EQ(permutations.get().size(),Combinatorics::factorial(s.size()));
    ASSERT_THAT(permutations.get()[0], ElementsAre(0));
}


TEST(ACombinatoricsTest, Permutations) {
    std::vector<int> s{0, 1, 2};

    Combinatorics::Permutations<int> permutations(s);

    ASSERT_EQ(permutations.get().size(),Combinatorics::factorial(s.size()));

    ASSERT_THAT(permutations.get()[0], ElementsAre(0, 1, 2));
    ASSERT_THAT(permutations.get()[1], ElementsAre(0, 2, 1));
    ASSERT_THAT(permutations.get()[2], ElementsAre(1, 0, 2));
    ASSERT_THAT(permutations.get()[3], ElementsAre(1, 2, 0));
    ASSERT_THAT(permutations.get()[4], ElementsAre(2, 0, 1));
    ASSERT_THAT(permutations.get()[5], ElementsAre(2, 1, 0));

}