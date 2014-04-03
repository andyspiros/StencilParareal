#pragma once

#include <cassert>
#include <vector>
#include <boost/static_assert.hpp>
#include "SharedInfrastructure.h"

/**
* @struct BlockConfiguration
* Structure holding block offset and size information
*/
struct BlockConfiguration
{
    int iBlockIndex;
    int iStart;
    int iEnd;
    int jBlockIndex;
    int jStart;
    int jEnd;
};

/**
* @class LaunchConfiguration
* Class storing the stencil launch information needed to apply a stencil with a certain boundary
*/
template<
    typename TBlockSize,
    typename TParameterIJBoundary,
    typename TLaunchConfigurationImpl>
class LaunchConfiguration 
{
    BOOST_STATIC_ASSERT(is_block_size<TBlockSize>::value);
    BOOST_STATIC_ASSERT(is_data_field_ij_boundary<TParameterIJBoundary>::value);
protected: // hide the constructors
    __CPU__
    LaunchConfiguration() {}
    __CPU__
    ~LaunchConfiguration() {}
    
    // method used by the actual launch configuration implementation to implement copy and assign
    __CPU__
    void Assign(const LaunchConfiguration& other) 
    {
        boundary_ = other.boundary_;
        blockConfigurations_ = other.blockConfigurations_;
    }

public:
    /**
    * Initialize the stencil launch configuration
    * @param calculationDomain calculation domain size
    * @param boundary boundary size
    */
    __CPU__
    void Init(const IJKSize& calculationDomain, const IJBoundary& boundary)
    {
        // copy boundary
        boundary_ = boundary;

        // make sure the boundary is smaller than or equal to the data field parameter boundaries
        assert(
            TParameterIJBoundary::IMinusOffset::value <= boundary.iMinusOffset() && 
            TParameterIJBoundary::IPlusOffset::value >= boundary.iPlusOffset() 
        );
        assert(
            TParameterIJBoundary::JMinusOffset::value <= boundary.jMinusOffset() && 
            TParameterIJBoundary::JPlusOffset::value >= boundary.jPlusOffset() 
        );

        // check that calculation domain and boundary are valid
        assert(calculationDomain.iSize() - boundary.iMinusOffset() + boundary.iPlusOffset() > 0);
        assert(calculationDomain.jSize() - boundary.jMinusOffset() + boundary.jPlusOffset() > 0); 
                
        // compute the block start index
        int iBlockStartIndex = 
            boundary_.iMinusOffset() < 0 ? 
            (boundary_.iMinusOffset() - TBlockSize::ISize::value + 1) / TBlockSize::ISize::value : 
            boundary_.iMinusOffset() / TBlockSize::ISize::value;
        int jBlockStartIndex = 
            boundary_.jMinusOffset() < 0 ? 
            (boundary_.jMinusOffset() - TBlockSize::JSize::value + 1) / TBlockSize::JSize::value : 
            boundary_.jMinusOffset() / TBlockSize::JSize::value;
        
        // compute the number of blocks
        int iNumberOfBlocks = 
            (calculationDomain.iSize() - iBlockStartIndex * TBlockSize::ISize::value + boundary_.iPlusOffset() + TBlockSize::ISize::value - 1) / 
            TBlockSize::ISize::value;
        int jNumberOfBlocks = 
            (calculationDomain.jSize() - jBlockStartIndex * TBlockSize::JSize::value + boundary_.jPlusOffset() + TBlockSize::JSize::value - 1) / 
            TBlockSize::JSize::value;
        
        // initialize block configuration vector
        blockConfigurations_.clear();
        // iterate over all blocks
        for(int jBlockIndex = jBlockStartIndex; jBlockIndex < jNumberOfBlocks + jBlockStartIndex; ++jBlockIndex)
        {
            for(int iBlockIndex = iBlockStartIndex; iBlockIndex < iNumberOfBlocks + iBlockStartIndex; ++iBlockIndex)
            {
                // compute the block start indexes which are usually 0
                int iStart = 0;
                int jStart = 0;
                
                // set block start values at minus boundary
                if(iBlockIndex * TBlockSize::ISize::value < boundary_.iMinusOffset())
                {
                    iStart = boundary.iMinusOffset() - iBlockIndex * TBlockSize::ISize::value;
                }
                if(jBlockIndex * TBlockSize::JSize::value < boundary_.jMinusOffset())
                {
                    jStart = boundary.jMinusOffset() - jBlockIndex * TBlockSize::JSize::value;
                }                
                
                // verify the start indexes
                assert(iStart >= 0 && iStart < TBlockSize::ISize::value);
                assert(jStart >= 0 && jStart < TBlockSize::JSize::value);

                // compute the block end indexes which are usually block size
                int iEnd = TBlockSize::ISize::value;
                int jEnd = TBlockSize::JSize::value;
            
                // set block start values at plus boundary
                if((iBlockIndex + 1) * TBlockSize::ISize::value > calculationDomain.iSize() + boundary.iPlusOffset())
                {
                    iEnd = calculationDomain.iSize() + boundary.iPlusOffset() - iBlockIndex * TBlockSize::ISize::value;
                }    
                if((jBlockIndex + 1) * TBlockSize::JSize::value > calculationDomain.jSize() + boundary.jPlusOffset())
                {
                    jEnd = calculationDomain.jSize() + boundary.jPlusOffset() - jBlockIndex * TBlockSize::JSize::value;
                }

                // verify the end indexes
                assert(iEnd > 0 && iEnd <= TBlockSize::ISize::value);
                assert(jEnd > 0 && jEnd <= TBlockSize::JSize::value);

                // initialize the block configuration
                BlockConfiguration blockConfiguration;
                blockConfiguration.iBlockIndex = iBlockIndex;
                blockConfiguration.iStart = iStart;
                blockConfiguration.iEnd = iEnd;
                blockConfiguration.jBlockIndex = jBlockIndex;
                blockConfiguration.jStart = jStart;
                blockConfiguration.jEnd = jEnd;
                 
                // store configuration in vector
                blockConfigurations_.push_back(blockConfiguration);
            }
        }
        
        // call back end specific init method
        static_cast<TLaunchConfigurationImpl*>(this)->InitImpl(blockConfigurations_);
    }

    /**
    * @return boundary of the stencil launch configuration
    */
    __CPU__
    const IJBoundary& boundary() const { return boundary_; }

    /**
    * @return vector of block configurations
    */ 
    __CPU__
    const std::vector<BlockConfiguration>& blockConfigurations() const { return blockConfigurations_; } 

private:
    IJBoundary boundary_;
    std::vector<BlockConfiguration> blockConfigurations_;
};

