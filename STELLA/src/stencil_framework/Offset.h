#pragma once

#include <boost/mpl/integral_c.hpp>

/**
* @struct Offset
* Type storing an access offset (used to access function and data parameters)
*/
template<int VI, int VJ, int VK>
struct Offset 
{
    typedef Offset<VI,VJ,VK> type;

    typedef boost::mpl::integral_c<int, VI> I;
    typedef boost::mpl::integral_c<int, VJ> J;
    typedef boost::mpl::integral_c<int, VK> K;
};

/**
* @struct add_access_offsets
* Meta function adding two access offsets
*/
template<
    typename TOffset1, 
    typename TOffset2>
struct add_offsets;

template<
    int VI1, 
    int VJ1, 
    int VK1, 
    int VI2, 
    int VJ2, 
    int VK2>
struct add_offsets<Offset<VI1, VJ1, VK1>, Offset<VI2, VJ2, VK2> >
{
    typedef Offset<VI1 + VI2, VJ1 + VJ2, VK1 + VK2> type; 
};

/**
* @struct scale_offset
* Meta function scaling an offset by a scaling factor
*/
template<
    typename TOffset,
    int VScalingFactor>
struct scale_offset;

template<
    int VI, 
    int VJ, 
    int VK, 
    int VScalingFactor>
struct scale_offset<Offset<VI, VJ, VK>, VScalingFactor>
{
    typedef Offset<VI * VScalingFactor, VJ * VScalingFactor, VK * VScalingFactor> type; 
};

// predefined offset and direction variables
#ifdef __CUDA_BACKEND__

typedef Offset<1, 0, 0> IDirection;
typedef Offset<0, 1, 0> JDirection;
__device__ static const IDirection idir;
__device__ static const JDirection jdir;

__device__ static const Offset<0, 0, 0> center;

__device__ static const Offset<1, 0, 0> iplus1;
__device__ static const Offset<2, 0, 0> iplus2;
__device__ static const Offset<3, 0, 0> iplus3;
__device__ static const Offset<-1, 0, 0> iminus1;
__device__ static const Offset<-2, 0, 0> iminus2;
__device__ static const Offset<-3, 0, 0> iminus3;

__device__ static const Offset<0, 1, 0> jplus1; 
__device__ static const Offset<0, 2, 0> jplus2; 
__device__ static const Offset<0, 3, 0> jplus3; 
__device__ static const Offset<0, -1, 0> jminus1;
__device__ static const Offset<0, -2, 0> jminus2;
__device__ static const Offset<0, -3, 0> jminus3;

__device__ static const Offset<0, 0, 1> kplus1;
__device__ static const Offset<0, 0, 2> kplus2;
__device__ static const Offset<0, 0, 3> kplus3;
__device__ static const Offset<0, 0, -1> kminus1;
__device__ static const Offset<0, 0, -2> kminus2;
__device__ static const Offset<0, 0, -3> kminus3;

#else

typedef Offset<1, 0, 0> IDirection;
typedef Offset<0, 1, 0> JDirection;
static const IDirection idir = IDirection();
static const JDirection jdir = JDirection();
       
static const Offset<0, 0, 0> center = Offset<0, 0, 0>();
       
static const Offset<1, 0, 0> iplus1 = Offset<1, 0, 0>();
static const Offset<2, 0, 0> iplus2 = Offset<2, 0, 0>();
static const Offset<3, 0, 0> iplus3 = Offset<3, 0, 0>();
static const Offset<-1, 0, 0> iminus1 = Offset<-1, 0, 0>();
static const Offset<-2, 0, 0> iminus2 = Offset<-2, 0, 0>();
static const Offset<-3, 0, 0> iminus3 = Offset<-3, 0, 0>();
       
static const Offset<0, 1, 0> jplus1 = Offset<0, 1, 0>(); 
static const Offset<0, 2, 0> jplus2 = Offset<0, 2, 0>(); 
static const Offset<0, 3, 0> jplus3 = Offset<0, 3, 0>(); 
static const Offset<0, -1, 0> jminus1 = Offset<0, -1, 0>();
static const Offset<0, -2, 0> jminus2 = Offset<0, -2, 0>();
static const Offset<0, -3, 0> jminus3 = Offset<0, -3, 0>();
       
static const Offset<0, 0, 1> kplus1 = Offset<0, 0, 1>();
static const Offset<0, 0, 2> kplus2 = Offset<0, 0, 2>();
static const Offset<0, 0, 3> kplus3 = Offset<0, 0, 3>();
static const Offset<0, 0, -1> kminus1 = Offset<0, 0, -1>();
static const Offset<0, 0, -2> kminus2 = Offset<0, 0, -2>();
static const Offset<0, 0, -3> kminus3 = Offset<0, 0, -3>();

#endif
  
