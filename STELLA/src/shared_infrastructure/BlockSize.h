#pragma once

#include <boost/static_assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>

/**
* @struct BlockSize
* Storing the size of a block
*/
template<
    int VISize, 
    int VJSize>
struct BlockSize
{
    BOOST_STATIC_ASSERT(VISize >= 0 && VJSize >= 0);

    typedef boost::mpl::integral_c<int, VISize> ISize;
    typedef boost::mpl::integral_c<int, VJSize> JSize;
    typedef boost::mpl::integral_c<int, VISize * VJSize> Size;
};

/**
* @struct is_block_size
* Meta function returning true if the parameter is a BlockSize
*/
template<typename T>
struct is_block_size : boost::mpl::false_ {};

template<int VISize, int VJSize>
struct is_block_size<BlockSize<VISize, VJSize> > : boost::mpl::true_ {};

  
