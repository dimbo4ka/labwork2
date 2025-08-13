#include <lib/number.h>
#include <gtest/gtest.h>
#include <bitset>
#include <cstring>
#include <tuple>


using TValue = std::pair<const char*, uint32_t>;

class ConvertingTestsSuite
    : public testing::TestWithParam<
        std::tuple<
            uint32_t,    // int value
            TValue,      // string value            
            bool         // equal
        >
    >
{
};

TEST_P(ConvertingTestsSuite, EqualTest) {
    uint239_t a = FromInt(std::get<0>(GetParam()), 0);
    uint239_t b = FromString(std::get<1>(GetParam()).first, std::get<1>(GetParam()).second);

    if(std::get<2>(GetParam()))
        ASSERT_EQ(a,b) << std::get<0>(GetParam()) << " == " << std::get<1>(GetParam()).first;
    else
        ASSERT_NE(a,b) << std::get<0>(GetParam()) << " != " << std::get<1>(GetParam()).first;
}


INSTANTIATE_TEST_SUITE_P(
    Group,
    ConvertingTestsSuite,
    testing::Values(
        // positive
        std::make_tuple(0,  TValue{"0", 0}, true),
        std::make_tuple(239, TValue{"239", 1}, true),
        std::make_tuple(255, TValue{"255", 2}, true),
        std::make_tuple(256, TValue{"256", 3}, true),
        std::make_tuple(10000, TValue{"10000", 4}, true),
        std::make_tuple(32767, TValue{"32767", 5}, true),
        std::make_tuple(32768, TValue{"32768", 6}, true),
        std::make_tuple(65535, TValue{"65535", 7}, true),
        std::make_tuple(2147483647, TValue{"2147483647", 8}, true),
        

        // negative
        std::make_tuple(0, TValue{"1", 0}, false),
        std::make_tuple(32768, TValue{"32769",  1}, false),
        std::make_tuple(255, TValue{"256", 2}, false),
        std::make_tuple(256, TValue{"255", 3}, false),
        std::make_tuple(31251, TValue{"31252", 4}, false),
        std::make_tuple(2147483647, TValue{"2147483648",  5}, false)
    )
);



class BitTestSuite 
    : public testing::TestWithParam<
        std::tuple<
            TValue,      // value
            const char*  // bit representation
        >
    >   
{
};


TEST_P(BitTestSuite, EqualTest) {
    constexpr size_t kBitsInByte = 8;
    auto toBitset = [](uint239_t& value) ->  std::bitset<35 * kBitsInByte> {
        std::bitset<35 * kBitsInByte> bits;
    
        for (std::size_t i = 0; i < 35; ++i) {
            for (std::size_t j = 0; j < kBitsInByte; ++j) {
                bits[i * kBitsInByte + j] = ( (value.data[34 - i] >> j) & 1 );
            }
        }

        return bits;
    };

    uint239_t a = FromString(std::get<0>(GetParam()).first, std::get<0>(GetParam()).second);
    std::bitset<35 * kBitsInByte> expected(std::get<1>(GetParam()));

    ASSERT_EQ(toBitset(a), expected);
}


INSTANTIATE_TEST_SUITE_P(
    Group,
    BitTestSuite,
    testing::Values(
        std::make_tuple(TValue{"1", 0},"0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"),
        std::make_tuple(TValue{"2024", 2024},"0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001001111111010000100000001000000010000000100000000000000010000000000000000000000000000000"),
        std::make_tuple(TValue{"1000000000000", 10},"0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101101000011010100010100110100010000000001000000000000000"),
        std::make_tuple(TValue{"777777777777777777", 77},"0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101001100101010011100001100001001001011101010000110000111000011100010000000000000000000000000000000010000000000000000000000010000000100000000000000010000000")
    )
);


class OperationTestsSuite
    : public testing::TestWithParam<
        std::tuple<
            TValue, // lhs
            TValue, // rhs
            TValue, // + result
            TValue, // - result
            TValue, // * result
            TValue // / result
        >
    > 
{
};


TEST_P(OperationTestsSuite, AddTest) {
    uint239_t a = FromString(std::get<0>(GetParam()).first, std::get<0>(GetParam()).second);
    uint239_t b = FromString(std::get<1>(GetParam()).first, std::get<1>(GetParam()).second);

    uint239_t result = a + b;
    uint239_t expected = FromString(std::get<2>(GetParam()).first, std::get<2>(GetParam()).second);

    ASSERT_EQ(result, expected);
}

TEST_P(OperationTestsSuite, SubstTest) {
    uint239_t a = FromString(std::get<0>(GetParam()).first, std::get<0>(GetParam()).second);
    uint239_t b = FromString(std::get<1>(GetParam()).first, std::get<1>(GetParam()).second);

    uint239_t result = a - b;
    uint239_t expected = FromString(std::get<3>(GetParam()).first, std::get<3>(GetParam()).second);

    ASSERT_EQ(result, expected);
}


TEST_P(OperationTestsSuite, MultTest) {
    uint239_t a = FromString(std::get<0>(GetParam()).first, std::get<0>(GetParam()).second);
    uint239_t b = FromString(std::get<1>(GetParam()).first, std::get<1>(GetParam()).second);

    uint239_t result = a * b;
    uint239_t expected = FromString(std::get<4>(GetParam()).first, std::get<4>(GetParam()).second);

    ASSERT_EQ(result, expected);
}

TEST_P(OperationTestsSuite, DivTest) {
    uint239_t a = FromString(std::get<0>(GetParam()).first, std::get<0>(GetParam()).second);
    uint239_t b = FromString(std::get<1>(GetParam()).first, std::get<1>(GetParam()).second);

    if(strcmp(std::get<1>(GetParam()).first, "0")){

        uint239_t result = a / b;
        uint239_t expected = FromString(std::get<5>(GetParam()).first, std::get<5>(GetParam()).second);

        ASSERT_EQ(result, expected);
    }
}


INSTANTIATE_TEST_SUITE_P(
    Group,
    OperationTestsSuite,
    testing::Values(
        std::make_tuple(TValue{"1", 0}, TValue{"1", 0}, TValue{"2", 0}, TValue{"0", 0},  TValue{"1", 0}, TValue{"1", 0}),
        std::make_tuple(TValue{"1", 2}, TValue{"1", 1}, TValue{"2", 3}, TValue{"0", 1},  TValue{"1", 3}, TValue{"1", 1}),
        std::make_tuple(TValue{"2024", 2024}, TValue{"8", 239}, TValue{"2032", 2263}, TValue{"2016", 1785},  TValue{"16192", 2263}, TValue{"253", 1785}),
        std::make_tuple(TValue{"876", 123}, TValue{"124", 48}, TValue{"1000", 171}, TValue{"752", 75},  TValue{"108624", 171}, TValue{"7", 75}),
        std::make_tuple(TValue{"99999999999999999999", 99}, TValue{"1", 0}, TValue{"100000000000000000000", 99}, TValue{"99999999999999999998", 99},  TValue{"99999999999999999999", 99}, TValue{"99999999999999999999", 99}),
        std::make_tuple(TValue{"1000", 1000}, TValue{"2", 999}, TValue{"1002", 1999}, TValue{"998", 1},  TValue{"2000", 1999}, TValue{"500", 1})
    )
);