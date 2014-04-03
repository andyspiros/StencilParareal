#include "gtest/gtest.h"

// define block sizes
const int iBlockSize = 3;
const int jBlockSize = 7;

// use openmp implementation to test block logic!
#include "LaunchConfigurationOpenMP.h"
typedef LaunchConfigurationOpenMP<
    BlockSize<iBlockSize, jBlockSize>, 
    DataFieldIJBoundary<-10, 10, -10, 10> 
> LaunchConfigurationImpl;

// function checking if an index is inside a block
bool check_in_block(const BlockConfiguration& blockConfiguration, const int i, const int j)
{
    bool result = false;
    if( i >= blockConfiguration.iBlockIndex * iBlockSize + blockConfiguration.iStart && 
        i < blockConfiguration.iBlockIndex * iBlockSize + blockConfiguration.iEnd &&
        j >= blockConfiguration.jBlockIndex * jBlockSize + blockConfiguration.jStart && 
        j < blockConfiguration.jBlockIndex * jBlockSize + blockConfiguration.jEnd       )
    {
        result = true;
    }
    return result;
}

// count occurrences of a position in the block configuration (we expect 0 or 1)
int count_block_occurrences(const std::vector<BlockConfiguration>& blockConfigurations, const int i, const int j)
{
    int count = 0;
    for(std::vector<BlockConfiguration>::const_iterator iter = blockConfigurations.begin(); iter != blockConfigurations.end(); ++iter)
    {
        if(check_in_block(*iter, i, j)) count++;
    }
    return count;
}

void configuration_test_round(int iDomain, int iMinusBoundary, int iPlusBoundary, int jDomain, int jMinusBoundary, int jPlusBoundary)
{
    IJKSize calculationDomain;
    IJBoundary boundary;

    calculationDomain.Init(iDomain, jDomain, 0);
    boundary.Init(iMinusBoundary, iPlusBoundary, jMinusBoundary, jPlusBoundary);

    LaunchConfigurationImpl launchConfig;
    launchConfig.Init(calculationDomain, boundary);

    // check the whole domain with a boundary 
    for(int i = iMinusBoundary - iBlockSize - 1; i < iDomain + iPlusBoundary + iBlockSize + 1; ++i)
    {
        for(int j = jMinusBoundary - jBlockSize - 1; j < iDomain + jPlusBoundary + jBlockSize + 1; ++j)
        {
            // check inside is covered
            if( i >= iMinusBoundary && i < iDomain + iPlusBoundary &&
                j >= jMinusBoundary && j < jDomain + jPlusBoundary    )
            {
                ASSERT_EQ(1, count_block_occurrences(launchConfig.blockConfigurations(), i, j));
            }
            // check boundary is not covered
            else
            {
                ASSERT_EQ(0, count_block_occurrences(launchConfig.blockConfigurations(), i, j));
            }
        }
    }
}

// test computation of block configuration vector
TEST(LaunchConfigurationUnittest, Init)
{
    // test standard case (calculation domain only multiple of block size)
    configuration_test_round(6, 0, 0, 7, 0, 0);    

    // inner boundary 
    configuration_test_round(5, 1, -3, 7, 2, -1);    

    // outer boundary
    configuration_test_round(5, -1, 3, 7, -7, 1);  

    // top left corner
    configuration_test_round(5, -1, -4, 7, -1, -6);  

    // top right corner
    configuration_test_round(5, 4, 2, 7, -1, -6);  

    // bottom left corner
    configuration_test_round(5, -1, -4, 7, 6, 8);  

    // bottom right corner
    configuration_test_round(5, 4, 2, 7, 6, 8);  

    // outside domain
    configuration_test_round(5, -5, -6, 7, -2, -8);  

    // small (1x1) inside domain 
    configuration_test_round(5, 3, -1, 7, 4, -2);  
}

// test copy
TEST(LaunchConfigurationUnittest, Copy)
{
    IJKSize calculationDomain;
    IJBoundary boundary;
    
    calculationDomain.Init(5, 7, 0);
    boundary.Init(-1, 3, -7, 1);
    
    LaunchConfigurationImpl launchConfig1, launchConfig2;
    launchConfig1.Init(calculationDomain, boundary);

    launchConfig2 = launchConfig1;

    // check boundary
    ASSERT_EQ(launchConfig1.boundary().iMinusOffset(), launchConfig2.boundary().iMinusOffset());
    ASSERT_EQ(launchConfig1.boundary().iPlusOffset(), launchConfig2.boundary().iPlusOffset());
    ASSERT_EQ(launchConfig1.boundary().jMinusOffset(), launchConfig2.boundary().jMinusOffset());
    ASSERT_EQ(launchConfig1.boundary().jPlusOffset(), launchConfig2.boundary().jPlusOffset());
    
    // check blocks
    ASSERT_EQ(launchConfig1.blockConfigurations().size(), launchConfig2.blockConfigurations().size());
    for(int i = 0; i < static_cast<int>(launchConfig1.blockConfigurations().size()); ++i)
    {
        ASSERT_EQ(launchConfig1.blockConfigurations()[i].iBlockIndex, launchConfig2.blockConfigurations()[i].iBlockIndex);
        ASSERT_EQ(launchConfig1.blockConfigurations()[i].iStart, launchConfig2.blockConfigurations()[i].iStart);
        ASSERT_EQ(launchConfig1.blockConfigurations()[i].iEnd, launchConfig2.blockConfigurations()[i].iEnd);
        ASSERT_EQ(launchConfig1.blockConfigurations()[i].jBlockIndex, launchConfig2.blockConfigurations()[i].jBlockIndex);
        ASSERT_EQ(launchConfig1.blockConfigurations()[i].jStart, launchConfig2.blockConfigurations()[i].jStart);
        ASSERT_EQ(launchConfig1.blockConfigurations()[i].jEnd, launchConfig2.blockConfigurations()[i].jEnd);
    }
}

