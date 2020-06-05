#include <gmock/gmock.h>

int main(int ac, char* av[])
{
    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
