#pragma once

#include <boost/config.hpp> 
#include <boost/static_assert.hpp>
#include <boost/type_traits/remove_pointer.hpp> 
#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/max.hpp>
#include "Enums.h"
#include "StencilStageEnvironment.h"
#include "StencilSweepDescriptor.h"
#include "StencilSweepGroupDescriptor.h"
#include "KLoopCascade.h"
#include "ContextCUDA.h"

/**
* @struct warps_per_row
* Meta function computing the number of warps per block row
*/
template<typename TBlockSize>
struct warps_per_row :
    boost::mpl::integral_c<int, (TBlockSize::ISize::value + cWarpSize - 1) / cWarpSize>  
{};

/**
* @struct warps_per_block
* Meta function computing the number of warps per block
*/
template<typename TBlockSize>
struct warps_per_block :
    boost::mpl::integral_c<int, TBlockSize::JSize::value * warps_per_row<TBlockSize>::value> 
{};

/**
* @struct warps_per_jminusboundary
* Meta function computing the number of warps per j minus boundary
*/
template<
    typename TBlockSize, 
    typename TIJRange>
struct warps_per_jminusboundary :
    boost::mpl::integral_c<int, -TIJRange::JMinusOffset::value * warps_per_row<TBlockSize>::value> 
{};

/**
* @struct warps_per_jplusboundary
* Meta function computing the number of warps per j plus boundary
*/
template<
    typename TBlockSize, 
    typename TIJRange>
struct warps_per_jplusboundary :
    boost::mpl::integral_c<int, TIJRange::JPlusOffset::value * warps_per_row<TBlockSize>::value> 
{};

/**
* @struct padded_boundary
* Meta function computing the boundary size padded to a power of 2 for an efficient split up of i boundary warps 
* (supports only boundaries up to size 4)
*/
template<int VBoundary>
struct padded_boundary :
    boost::mpl::integral_c<int, VBoundary <= 1 ? 1 : (VBoundary <= 2 ? 2 : 4)> 
{
    BOOST_STATIC_ASSERT(VBoundary <= 4);
};

/**
* @struct rows_per_iboundary_warp
* Meta function computing the number of rows per i boundary warp
*/
template<int VBoundary>
struct rows_per_iboundary_warp :
    boost::mpl::integral_c<int, cWarpSize / padded_boundary<VBoundary>::value> 
{
    BOOST_STATIC_ASSERT(cWarpSize % padded_boundary<VBoundary>::value == 0);
};

/**
* @struct iboundaries_jsize 
* Meta function computing the j size of the i boundary 
*/
template<
    typename TBlockSize, 
    typename TIJRange>
struct iboundaries_jsize :
    boost::mpl::integral_c<int, TBlockSize::JSize::value + 
        (
            TIJRange::CornerPolicy::value == cIndented ? 
            0 : 
            TIJRange::JPlusOffset::value - TIJRange::JMinusOffset::value
        )
    > 
{};

/**
* @struct iboundaries_jstart 
* Meta function computing the j start offset of the i boundaries
*/
template<
    typename TBlockSize, 
    typename TIJRange>
struct iboundaries_jstart :
    boost::mpl::integral_c<int, TIJRange::CornerPolicy::value == cIndented ? 0 : TIJRange::JMinusOffset::value> 
{};

/**
* @struct warps_per_iminusboundary
* Meta function computing the number of warps per i minus boundary
*/
template<
    typename TBlockSize, 
    typename TIJRange>
struct warps_per_iminusboundary :
    boost::mpl::integral_c<int, 
        TIJRange::IMinusOffset::value == 0 ? 
        0 :
        (iboundaries_jsize<TBlockSize, TIJRange>::value + rows_per_iboundary_warp<TIJRange::IMinusOffset::value>::value - 1) / 
        rows_per_iboundary_warp<TIJRange::IMinusOffset::value>::value
    >
{};

/**
* @struct warps_per_iplusboundary
* Meta function computing the number of warps per i plus boundary
*/
template<
    typename TBlockSize, 
    typename TIJRange>
struct warps_per_iplusboundary :
    boost::mpl::integral_c<int, 
        TIJRange::IPlusOffset::value == 0 ? 
        0 :
        (iboundaries_jsize<TBlockSize, TIJRange>::value + rows_per_iboundary_warp<TIJRange::IPlusOffset::value>::value - 1) /
        rows_per_iboundary_warp<TIJRange::IPlusOffset::value>::value
    >
{};

/**
* @struct warps_per_sweep
* Meta function computing the number of warps needed for a given sweep descriptor
*/
template<
    typename TBlockSize, 
    typename TStencilSweepDescriptor>
struct warps_per_sweep 
{
    typedef typename stencil_sweep_descriptor_maximum_ij_range<TStencilSweepDescriptor>::type IJRange;
    typedef boost::mpl::integral_c<int, 
        (
            warps_per_block<TBlockSize>::value +
            warps_per_jminusboundary<TBlockSize, IJRange>::value +
            warps_per_jplusboundary<TBlockSize, IJRange>::value +
            warps_per_iminusboundary<TBlockSize, IJRange>::value +
            warps_per_iplusboundary<TBlockSize, IJRange>::value 
        )
    > type;
};

/**
* @struct warps_per_sweep_group
* Meta function computing the number of warps needed for a given sweep group descriptor
*/
template<
    typename TBlockSize, 
    typename TStencilSweepGroupDescriptor>
struct warps_per_sweep_group
{
    // extract the stencil sweep descriptors of the group
    typedef typename stencil_sweep_group_descriptor_sweeps<
        TStencilSweepGroupDescriptor
    >::type StencilSweepDescriptors;

    // find the maximal number of warps needed for the different sweeps
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::fold<
                StencilSweepDescriptors,
                boost::mpl::integral_c<int,0>,
                boost::mpl::max<
                    boost::mpl::_1, 
                    warps_per_sweep<TBlockSize, boost::mpl::_2>
                >
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct warps_per_stencil
* Meta function computing the number of warps needed for a given stencil
*/
template<
    typename TBlockSize, 
    typename TStencilSweepGroupDescriptors>
struct warps_per_stencil 
{
    // find the maximal number of warps needed for the different groups
    BOOST_STATIC_CONSTANT(int, value = 
        (
            boost::mpl::fold<
                TStencilSweepGroupDescriptors,
                boost::mpl::integral_c<int,0>,
                boost::mpl::max<
                    boost::mpl::_1, 
                    warps_per_sweep_group<TBlockSize, boost::mpl::_2>
                >
            >::type::value
        )
    );
    typedef boost::mpl::integral_c<int, int(value)> type;
};

/**
* @struct BlockLoopFunctorCUDA
* Functor implementing the CUDA loop which applies a sweep to a block. The idea is that a certain amount 
* of the warps in a block update the inside of the block plus the j boundary, while the remaining warps are
* used to update i minus respectively i plus boundary. Note that there might be warps doing nothing
* at all. An i boundary warp handles multiple rows of the i boundary.
*/
struct BlockLoopFunctorCUDA
{
    template<
        typename TSharedDataPointer, 
        typename TStencilSweepDescriptor>
    __ACC__
    static void Do(typename parameter_type<TSharedDataPointer>::type pSharedData)
    {        
        // define the loop offsets  
        typedef typename stencil_sweep_descriptor_maximum_ij_range<TStencilSweepDescriptor>::type IJRange;
        
        // define the context type
        typedef typename boost::remove_pointer<TSharedDataPointer>::type SharedData;
        typedef ContextCUDA<SharedData, TStencilSweepDescriptor> Context; 
        typedef typename Context::BlockSize BlockSize;

        // define the k loop cascade
        typedef KLoopCascade<Context, TStencilSweepDescriptor> KLoop; 

        // make sure the i block size is a multiple of the warp size
        BOOST_STATIC_ASSERT(BlockSize::ISize::value % cWarpSize == 0);

        // compute the limits used to assign the warps to a certain region of the block / boundary
        const int blockAndJBoundaryLimit = 
            warps_per_block<BlockSize>::value +
            warps_per_jminusboundary<BlockSize, IJRange>::value + 
            warps_per_jplusboundary<BlockSize, IJRange>::value;

        const int iMinusBoundaryLimit = 
            blockAndJBoundaryLimit +
            warps_per_iminusboundary<BlockSize, IJRange>::value;

        const int iPlusBoundaryLimit = 
            iMinusBoundaryLimit +
            warps_per_iplusboundary<BlockSize, IJRange>::value;

        // setup the context
        Context context; 
        context.Init(pSharedData);

        // handle block inside and j boundary
        if(static_cast<int>(threadIdx.y) < blockAndJBoundaryLimit)
        {
            // line up the warps row by row starting from the bottom of the block / j boundary
            int jOffset = static_cast<int>(threadIdx.y) / warps_per_row<BlockSize>::value; // warp index divided by the number of warps per row
            const int iOffset =
                static_cast<int>(threadIdx.x) + // thread index inside the warp
                (static_cast<int>(threadIdx.y) - jOffset * warps_per_row<BlockSize>::value) * cWarpSize; // number of preceding warps in the row times number of thread per warp
            jOffset = jOffset + IJRange::JMinusOffset::value; // add the block plus j boundary start offset

            // setup the context position
            KLoop::RestoreAndAdvance(context, iOffset, jOffset);
        }
        // handle the i minus boundary
        else if(static_cast<int>(threadIdx.y) < iMinusBoundaryLimit)
        {
            // define a i minus boundary warp index starting from zero
            const int warpIndex = static_cast<int>(threadIdx.y) - blockAndJBoundaryLimit;

            // define the number of rows updated by a single wrap and the padded width of the warp
            const int boundaryRowsPerWarp = rows_per_iboundary_warp<-IJRange::IMinusOffset::value>::value;
            const int paddedBoundary = padded_boundary<-IJRange::IMinusOffset::value>::value;

            // split up the warps into multiple segments handling a single row of the i boundary (pad the row segment to a power of 2)
            const int iOffset = 
                -paddedBoundary + // start offset of the i minus boundary row segments 
                (static_cast<int>(threadIdx.x) % paddedBoundary); // row segment offset
            const int jOffset = 
                warpIndex * boundaryRowsPerWarp + // number of preceding warps times rows updated per warp
                static_cast<int>(threadIdx.x) / paddedBoundary +  // number of preceding row segments
                iboundaries_jstart<BlockSize, IJRange>::value; // start offset of the i boundary

            // setup the context position
            KLoop::RestoreAndAdvance(context, iOffset, jOffset);
        }
        // handle the i plus boundary
        else if(static_cast<int>(threadIdx.y) < iPlusBoundaryLimit)
        {
            // define a i plus boundary warp index starting from zero
            const int warpIndex = static_cast<int>(threadIdx.y) - iMinusBoundaryLimit;

            // define the number of rows updated by a single wrap and the padded width of the warp
            const int boundaryRowsPerWarp = rows_per_iboundary_warp<IJRange::IPlusOffset::value>::value;
            const int paddedBoundary = padded_boundary<IJRange::IPlusOffset::value>::value;
             
            // split up the warps into multiple row segments handling a single row of the i boundary (pad the row segment to a power of 2)
            const int iOffset = 
                static_cast<int>(threadIdx.x) % paddedBoundary + // row segment offset
                BlockSize::ISize::value;  // start offset of the i plus boundary row segments 
            const int jOffset = 
                warpIndex * boundaryRowsPerWarp + // number of preceding warps times rows updated per warp
                static_cast<int>(threadIdx.x) / paddedBoundary + // number of preceding row segments
                iboundaries_jstart<BlockSize, IJRange>::value; // start offset of the i boundary

            // setup the context position
            KLoop::RestoreAndAdvance(context, iOffset, jOffset);
        }
        else
        {
            // setup the context position outside of the IJRange 
            // (note that the stencil sweep functor will pause these threads)
            KLoop::RestoreAndAdvance(context, IJRange::IMinusOffset::value - 1, IJRange::JMinusOffset::value - 1);
        }

        // loop over k
        KLoop::Do(context);  

        // synchronize threads at the end of the k loop
        __syncthreads();
    }
};






