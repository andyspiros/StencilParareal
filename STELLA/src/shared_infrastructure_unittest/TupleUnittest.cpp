#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/vector/vector10_c.hpp>
#include <boost/mpl/vector/vector10.hpp>
#include <string>
#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

TEST(TupleUnittest, CreationAndAccess)
{
    std::string text = "hallo";
    typedef TupleElements<
        boost::mpl::vector3_c<int,3,5,7>, 
        boost::mpl::vector3<int,double,std::string> 
    > Elements;
    Tuple<Elements> tuple;

    tuple.Init(20, 13.5, text);

    ASSERT_EQ(20, (tuple(static_cast<boost::mpl::integral_c<int,3>*>(0))));
    ASSERT_EQ(13.5, (tuple(static_cast<boost::mpl::integral_c<int,5>*>(0))));
    ASSERT_EQ("hallo", (tuple(static_cast<boost::mpl::integral_c<int,7>*>(0))));
}

TEST(TupleUnittest, CopyConstructor)
{
    typedef Tuple<
        TupleElements<
            boost::mpl::vector2_c<int,7,2>, 
            boost::mpl::vector2<int,double>
        >
    > TupleType;

    TupleType orig;
    orig.Init(20, 13.5);

    TupleType copy = orig;

    ASSERT_EQ(copy(static_cast<boost::mpl::integral_c<int,7>*>(0)), orig(static_cast<boost::mpl::integral_c<int,7>*>(0)));
    ASSERT_EQ(copy(static_cast<boost::mpl::integral_c<int,2>*>(0)), orig(static_cast<boost::mpl::integral_c<int,2>*>(0)));
}

struct ConcatenateStrings
{
    template<typename T>
    static void Do(T&, parameter_type<std::string>::type) {}

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

TEST(TupleUnittest, ModifyTuple)
{
    std::string text = "hallo";
    typedef Tuple<
        TupleElements<
            boost::mpl::vector2_c<int,0,1>, 
            boost::mpl::vector2<int,std::string>
        >
    > TupleType;
    TupleType tuple;
    tuple.Init(20, text);
    
    std::string world = " world";
    modify_tuple<
        ConcatenateStrings,
        boost::mpl::vector2_c<int, 0, 1>,
        std::string
    >(tuple, world);

    ASSERT_EQ(20, ( tuple(static_cast<boost::mpl::integral_c<int, 0>*>(0)) ) );
    ASSERT_EQ("hallo world", ( tuple(static_cast<boost::mpl::integral_c<int, 1>*>(0)) ) );

    modify_tuple<
        ConcatenateStrings,
        boost::mpl::vector1_c<int, 0>,
        std::string
    >(tuple, world);

    // assert still just 
    ASSERT_EQ(20, ( tuple(static_cast<boost::mpl::integral_c<int, 0>*>(0)) ) );
    ASSERT_EQ("hallo world", ( tuple(static_cast<boost::mpl::integral_c<int, 1>*>(0)) ) );
    
    modify_tuple<
        AddExclamationMark,
        boost::mpl::vector1_c<int, 1>
    >(tuple);

    // now plus !
    ASSERT_EQ(20, ( tuple(static_cast<boost::mpl::integral_c<int, 0>*>(0)) ) );
    ASSERT_EQ("hallo world!", ( tuple(static_cast<boost::mpl::integral_c<int, 1>*>(0)) ) );
}

struct StringToStringLength
{
    static void Do(std::string& text, int& len)
    {
        len = static_cast<int>(text.length());
    }
};

struct StringToTotalStringLength
{
    static void Do(std::string& text, int& len, parameter_type<int&>::type sum)
    {
        len = static_cast<int>(text.length()) + sum;
        sum = len;
    }
};

TEST(TupleUnittest, Modify2Tuples)
{
    std::string text1 = "hallo1";
    std::string text2 = "hallo12";
    std::string text3 = "hallo123";

    typedef Tuple<
        TupleElements<
            boost::mpl::vector3_c<int,0,1,2>, 
            boost::mpl::vector3<std::string,std::string,std::string> 
        >
    > TupleType1;
    TupleType1 tuple1;
    tuple1.Init(text1, text2, text3);

    typedef Tuple<
        TupleElements<
            boost::mpl::vector3_c<int,0,1,2>, 
            boost::mpl::vector3<int,int,int> 
        >
    > TupleType2;
    TupleType2 tuple2;
    tuple2.Init(0, 0, 0);

    // check modify 2 tuples without parameter works
    modify_2_tuples<
        StringToStringLength,
        boost::mpl::vector2_c<int, 0, 2>
    >(tuple1, tuple2); 

    ASSERT_EQ(6, (tuple2(static_cast<boost::mpl::integral_c<int, 0>*>(0))));
    ASSERT_EQ(0, (tuple2(static_cast<boost::mpl::integral_c<int, 1>*>(0))));
    ASSERT_EQ(8, (tuple2(static_cast<boost::mpl::integral_c<int, 2>*>(0))));
    
    // check modify 2 tuples with parameter works
    int sum = 0;
    modify_2_tuples<
        StringToTotalStringLength,
        boost::mpl::vector2_c<int, 0, 2>,
        int&
    >(tuple1, tuple2, sum); 

    ASSERT_EQ(6, (tuple2(static_cast<boost::mpl::integral_c<int, 0>*>(0))));
    ASSERT_EQ(0, (tuple2(static_cast<boost::mpl::integral_c<int, 1>*>(0))));
    ASSERT_EQ(14, (tuple2(static_cast<boost::mpl::integral_c<int, 2>*>(0))));
    ASSERT_EQ(14, sum);
}


