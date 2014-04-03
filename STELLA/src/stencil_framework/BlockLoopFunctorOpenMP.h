#pragma once

#include <boost/mpl/integral_c.hpp>
#include "Enums.h"
#include "StencilStageEnvironment.h"
#include "StencilSweepDescriptor.h"
#include "KLoopCascade.h"

/**
* @struct BlockLoopFunctorOpenMP
* Functor implementing the OpenMP loop which applies a sweep to a block
*/
struct BlockLoopFunctorOpenMP 
{
    template<
        typename TContext, 
        typename TStencilSweepDescriptor>
    static void Do(typename parameter_type<TContext>::type context)
    {
        // choose the right implementation depending on the corner strategy
        doImpl<
            TContext,
            TStencilSweepDescriptor
        >(context, static_cast<typename stencil_sweep_descriptor_maximum_ij_range<TStencilSweepDescriptor>::type::CornerPolicy*>(0));
    }

private:
    // implement loop over column including the corners
    template<
        typename TContext, 
        typename TStencilSweepDescriptor>
    static void doImpl(TContext& context, boost::mpl::integral_c<CornerPolicy, cComplete>*)
    {
        // define the block size and the loop offsets  
        typedef typename TContext::BlockSize BlockSize;

        // calculate ij-loop boundaries
        typedef typename stencil_sweep_descriptor_maximum_ij_range<TStencilSweepDescriptor>::type IJRange;
        
        // define the k loop cascade
        typedef KLoopCascade<TContext, TStencilSweepDescriptor> KLoop; 
        
        // calculate ij loop boundaries
        const int iStart = context.iBlockStart() + IJRange::IMinusOffset::value;
        const int iEnd = context.iBlockEnd() + IJRange::IPlusOffset::value;
        const int jStart = context.jBlockStart() + IJRange::JMinusOffset::value;
        const int jEnd = context.jBlockEnd() + IJRange::JPlusOffset::value;

        // loop over the block
        for(int j = jStart; j < jEnd; ++j)
        {
            // set i start value
            int i = iStart;
            
            // restore the context position
            KLoop::RestoreAndAdvance(context, i, j);
            
            // loop over i
            while(true)
            {
                // loop over k
                KLoop::Do(context);

                // check if i loop finished if not goto the next i-column
                if(++i < iEnd) KLoop::template ColumnAdvance<1,0>(context);
                else break;
            }
        }
    }

    // implement loop over a column excluding the corners
    template<
        typename TContext, 
        typename TStencilSweepDescriptor>
    static void doImpl(TContext& context, boost::mpl::integral_c<CornerPolicy, cIndented>*)
    {
        // define the block size and the loop offsets  
        typedef typename TContext::BlockSize BlockSize;

        // calculate ij-loop boundaries
        typedef typename stencil_sweep_descriptor_maximum_ij_range<TStencilSweepDescriptor>::type IJRange;

        // define the k loop cascade
        typedef KLoopCascade<TContext, TStencilSweepDescriptor> KLoop;

        // calculate j loop boundaries
        const int jStart = context.jBlockStart() + IJRange::JMinusOffset::value;
        const int jEnd = context.jBlockEnd() + IJRange::JPlusOffset::value;

        // loop over the block
        for(int j = jStart; j < jEnd; ++j)
        {
            // compute if we are in boundary
            const bool inBoundary = (j < 0) || (j >= BlockSize::JSize::value);
            
            // compute i start and i end value
            const int iEnd = context.iBlockEnd() + (inBoundary ? 0 : IJRange::IPlusOffset::value);
            int i = context.iBlockStart() + (inBoundary ? 0 : IJRange::IMinusOffset::value);
            
            // restore the context position
            KLoop::RestoreAndAdvance(context, i, j);

            // loop over i
            while(true)
            {
                // loop over k
                KLoop::Do(context);

                // check if i loop finished if not goto the next i-column
                if(++i < iEnd) KLoop::template ColumnAdvance<1,0>(context);
                else break;
            }
        }
    }
};






