#include "gtest/gtest.h"
#include "SharedInfrastructure.h"

TEST(UsageMeterUnittest, All)
{
    const int numberOfIterations = 5;
    UsageMeter applyMeter;
    UsageMeter bytesMeter;

    applyMeter.Init("applyMeter", "counter");
    bytesMeter.Init("bandwidthMeter", "bytes");

    // test increment reset
    applyMeter.Add(2.3);
    ASSERT_EQ(applyMeter.counter(), 2.3);
    applyMeter.Reset();
    ASSERT_EQ(applyMeter.counter(), 0.0);

    for (int i=0; i<numberOfIterations; ++i)
    {
        applyMeter.Add(1.0);
        bytesMeter.Add(100.0);
    }
    ASSERT_EQ(applyMeter.counter(), 1.0*numberOfIterations);
    ASSERT_EQ(bytesMeter.counter(), 100.0*numberOfIterations);
}
