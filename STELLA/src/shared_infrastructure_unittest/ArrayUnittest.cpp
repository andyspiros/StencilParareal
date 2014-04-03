#include <string>
#include <boost/mpl/integral_c.hpp>
#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

TEST(ArrayUnittest, CreationAndAccess)
{
    Array<double, 2> doubleArray;
    Array<std::string, 3> stringArray;

    doubleArray.At(static_cast<boost::mpl::integral_c<int,0>*>(0)) = 2.3;
    doubleArray.At(static_cast<boost::mpl::integral_c<int,1>*>(0)) = 4.5;

    stringArray.At(static_cast<boost::mpl::integral_c<int,0>*>(0)) = "hallo";
    stringArray.At(static_cast<boost::mpl::integral_c<int,1>*>(0)) = " ";
    stringArray.At(static_cast<boost::mpl::integral_c<int,2>*>(0)) = "world";

    ASSERT_EQ(2.3, (doubleArray.At(static_cast<boost::mpl::integral_c<int,0>*>(0))));
    ASSERT_EQ(4.5, (doubleArray.At(static_cast<boost::mpl::integral_c<int,1>*>(0))));
    ASSERT_EQ("hallo", (stringArray.At(static_cast<boost::mpl::integral_c<int,0>*>(0))));
    ASSERT_EQ(" ", (stringArray.At(static_cast<boost::mpl::integral_c<int,1>*>(0))));
    ASSERT_EQ("world", (stringArray.At(static_cast<boost::mpl::integral_c<int,2>*>(0))));
}

TEST(ArrayUnittest, CopyConstructor)
{
    Array<double, 2> orig;
    Array<double, 2> copy;

    orig.At(static_cast<boost::mpl::integral_c<int,0>*>(0)) = 2.3;
    orig.At(static_cast<boost::mpl::integral_c<int,1>*>(0)) = 4.5;

    copy = orig;

    ASSERT_EQ(copy.At(static_cast<boost::mpl::integral_c<int,0>*>(0)), orig.At(static_cast<boost::mpl::integral_c<int,0>*>(0)));
    ASSERT_EQ(copy.At(static_cast<boost::mpl::integral_c<int,1>*>(0)), orig.At(static_cast<boost::mpl::integral_c<int,1>*>(0)));
}

struct ConcatenateStrings
{
    static void Do(std::string& elem, parameter_type<std::string>::type param)
    {
        elem += param;
    }
};

struct AddExclamationMark
{
    static void Do(std::string& elem)
    {
        elem += "!";
    }
};

TEST(ArrayUnittest, ModifyArray)
{
    Array<std::string, 2> array;
    array.At(static_cast<boost::mpl::integral_c<int, 0>*>(0)) = "hallo"; 
    array.At(static_cast<boost::mpl::integral_c<int, 1>*>(0)) = "ciao";

    std::string world = " world";
    modify_array<
        ConcatenateStrings,
        boost::mpl::vector2_c<int, 0, 1>,
        std::string
    >(array, world);

    ASSERT_EQ("hallo world", array.At(static_cast<boost::mpl::integral_c<int, 0>*>(0)));
    ASSERT_EQ("ciao world", array.At(static_cast<boost::mpl::integral_c<int, 1>*>(0)));

    modify_array<
        AddExclamationMark,
        boost::mpl::vector1_c<int, 1>
    >(array);

    // check that only ciao world was extended with an exclamation mark
    ASSERT_EQ("hallo world", array.At(static_cast<boost::mpl::integral_c<int, 0>*>(0)));
    ASSERT_EQ("ciao world!", array.At(static_cast<boost::mpl::integral_c<int, 1>*>(0)));
}
