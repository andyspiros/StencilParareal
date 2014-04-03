#include <boost/mpl/vector/vector10_c.hpp>
#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

struct SumData
{
    int sum;
    int terms;
};

struct MultiplyData
{
    int product;
};

struct SumFunctor
{
    template<
        typename TSumData,
        typename T>
    static void Do(typename parameter_type<TSumData>::type data)
    {
        data.sum += T::value;
        data.terms++;
    }
};

struct SumAndMultiplyFunctor
{
    template<
        typename TSumData,
        typename TMultiplyData,
        typename T>
    static void Do(typename parameter_type<TSumData>::type data1, typename parameter_type<TMultiplyData>::type data2)
    {
        data1.sum += T::value;
        data1.terms++;
        data2.product *= T::value;
    }
};

TEST(ApplyToAllUnittest, SingleDataParameter)
{
    SumData data;
    data.sum = 0;
    data.terms = 0;
    apply_to_all<SumFunctor, boost::mpl::vector3_c<int,1,2,3>, SumData>(data);
    
    // 1 + 2 + 3 = 6
    ASSERT_EQ(6, data.sum);
    ASSERT_EQ(3, data.terms);
    
    data.sum = 0;
    data.terms = 0;
    apply_to_all<SumFunctor, boost::mpl::vector4_c<int,3,6,8,15>, SumData>(data);

    // 3 + 6 + 8 + 15 = 32
    ASSERT_EQ(32, data.sum);
    ASSERT_EQ(4, data.terms);
}

TEST(ApplyToAllUnittest, TwoDataParameters)
{
    SumData data1;
    MultiplyData data2; 
    data1.sum = 0;
    data1.terms = 0;
    data2.product = 1;
    apply_to_all<SumAndMultiplyFunctor, boost::mpl::vector3_c<int,1,2,3>, SumData, MultiplyData>(data1, data2);

    // 1 + 2 + 3 = 6
    // 1 * 2 * 3 = 6
    ASSERT_EQ(6, data1.sum);
    ASSERT_EQ(3, data1.terms);
    ASSERT_EQ(6, data2.product);

    data1.sum = 0;
    data1.terms = 0;
    data2.product = 1;
    apply_to_all<SumAndMultiplyFunctor, boost::mpl::vector4_c<int,3,6,8,15>, SumData, MultiplyData>(data1, data2);

    // 3 + 6 + 8 + 15 = 32
    // 3 * 6 * 8 * 15 = 2160
    ASSERT_EQ(32, data1.sum);
    ASSERT_EQ(4, data1.terms);
    ASSERT_EQ(2160, data2.product);
}

