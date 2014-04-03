#pragma once

#include <cassert>
#include <string>
#include <boost/mpl/assert.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/transform.hpp>
#include "SharedInfrastructure.h"
#include "StencilSweepGroupFunctorOpenMP.h"
#include "ContextOpenMP.h"
#include "ColumnBufferOpenMP.h"
#include "LaunchConfigurationOpenMP.h"
#include "LaunchConfigurationManager.h"
#include "Stencil.h"
#include "StencilOpenMPFunctors.h"

#ifdef __OPENMP_BACKEND__
#include <omp.h>
#endif

/**
* @class StencilWrapperOpenMP
* OpenMP implementation of the stencil interface
*/
template<
    typename TStencilConfiguration,
    typename TParameterIJBoundary,
    typename TParameterTupleElements,
    typename TTemporaryFields,
    typename TStencilSweepGroupDescriptors>
class StencilOpenMP : public StencilInterface // implements the stencil interface
{
    DISALLOW_COPY_AND_ASSIGN(StencilOpenMP);
public:
    // define the block size
    typedef typename TStencilConfiguration::BlockSize BlockSize;
    typedef LaunchConfigurationOpenMP<BlockSize, TParameterIJBoundary> LaunchConfiguration;

    // define the parameter tuple information
    typedef typename TParameterTupleElements::ElementIndexes ParameterElementIndexes;
    typedef typename TParameterTupleElements::ElementTypes ParameterElementTypes;
    typedef TParameterTupleElements ParameterTupleElements;
    typedef Tuple<ParameterTupleElements> ParameterTupleType;

    // calculate the indexes of all data field parameters
    typedef typename data_field_parameters<ParameterElementTypes>::type DataFieldParameterIndexes;

    // define the temporary field tuple information
    typedef typename boost::mpl::transform<
        TTemporaryFields,
        temporary_field_index<boost::mpl::_>
    >::type TemporaryFieldElementIndexes;
    typedef typename boost::mpl::transform<
        TTemporaryFields,
        create_temporary_field<boost::mpl::_, BlockSize, TParameterIJBoundary, TStencilSweepGroupDescriptors>
    >::type TemporaryFieldElementTypes;
    typedef TupleElements<TemporaryFieldElementIndexes, TemporaryFieldElementTypes> TemporaryFieldTupleElements;
    typedef Tuple<TemporaryFieldTupleElements> TemporaryFieldTupleType;

    // define the apply tuple
    typedef typename merge_tuple_elements<
        ParameterTupleElements,
        TemporaryFieldTupleElements,
        parameter_to_apply<boost::mpl::_>,
        temporary_field_to_apply<boost::mpl::_>
    >::type ApplyTupleElements;

    // define the context descriptor type
    typedef ContextDescriptor<
        TStencilConfiguration, 
        ApplyTupleElements
    > ContextDescriptorType;
    typedef ContextOpenMP<ContextDescriptorType> Context;

    StencilOpenMP() {}
    virtual ~StencilOpenMP() {}

    /**
    * Initialize the wrapper
    * @param calculationDomain size of the calculation domain
    * @param parameterTuple tuple holding the stencil parameters
    */
    void Init(const IJKSize& calculationDomain, Tuple<ParameterTupleElements>& parameterTuple)
    {
        // setup the launch configuration manager
        IJBoundary defaultBoundary;
        defaultBoundary.Init(0, 0, 0, 0);
        launchConfigurationManager_.Init(calculationDomain, defaultBoundary);

        // initialize the parameter tuple
        parameterTuple_ = parameterTuple;

#ifdef __OPENMP_BACKEND__ 
        // setup thread private contexts
        threadPrivateContexts_.resize(omp_get_max_threads());
        threadPrivateTemporaryFieldTuples_.resize(omp_get_max_threads());
        #pragma omp parallel
        {
            int index = omp_get_thread_num();
            assert(static_cast<int>(threadPrivateContexts_.size()) > index);
            assert(static_cast<int>(threadPrivateTemporaryFieldTuples_.size()) > index);
#else
        // if not compiling with OpenMP setup a single thread private context
        threadPrivateContexts_.resize(1);
        threadPrivateTemporaryFieldTuples_.resize(1);
        {
            int index = 0;
#endif           
            // initialize the temporary field tuple
            modify_tuple< 
                TemporaryFieldInitFunctor,
                TemporaryFieldElementIndexes,
                const IJKSize
            >(threadPrivateTemporaryFieldTuples_[index], calculationDomain);

            // initialize the apply tuple
            threadPrivateContexts_[index].Init(calculationDomain.kSize());
            modify_2_tuples<
                ParameterToApplyFunctor,
                ParameterElementIndexes
            >(parameterTuple_, threadPrivateContexts_[index].dataTuple());

            modify_2_tuples<
                TemporaryFieldToApplyFunctor,
                TemporaryFieldElementIndexes
            >(threadPrivateTemporaryFieldTuples_[index], threadPrivateContexts_[index].dataTuple());
        }
    }

    // stencil interface implementation

    /**
    * Method applying the stencil using the default boundary
    */
    virtual void ApplyImpl()
    {
        // run the stencil using the default configuration
        applyImpl(launchConfigurationManager_.defaultConfiguration());
    }

    /**
    * Method applying the stencil using a variable boundary
    * @param boundary apply the stencil with a boundary around the calculation domain
    */
    virtual void ApplyImpl(const IJBoundary& boundary)
    {
        // find a matching launch configuration
        const LaunchConfiguration* pLaunchConfiguration = launchConfigurationManager_.FindConfiguration(boundary);

        // if a config was found apply it 
        if(pLaunchConfiguration) 
        { 
            applyImpl(*pLaunchConfiguration);
        }
        // if no config was found add a new one and apply it
        else
        {
            applyImpl(launchConfigurationManager_.AddConfiguration(boundary));
        }
    }

    /**
    * Method printing the bounding box for all stencil parameters
    */
    virtual std::string PrintBoundingBoxInfoImpl()
    {
#ifndef NDEBUG
        // update the pointers
        assert(threadPrivateContexts_.size() > 0);
        std::ostringstream boundingBoxInfo;
        modify_2_tuples<
            PrintBoundingBoxFunctor,
            DataFieldParameterIndexes,
            std::ostringstream
        >(parameterTuple_, threadPrivateContexts_[0].dataTuple(), boundingBoxInfo);

        return boundingBoxInfo.str();
#else
        return "PrintBoundingBoxInfo not supported";
#endif
    }

private:
    // method implementing apply given a launch configuration
    void applyImpl(const LaunchConfiguration& launchConfiguration)
    {
        const int numberOfBlocks = static_cast<int>(launchConfiguration.blockConfigurations().size());
        if(numberOfBlocks > 0)
        {    
#ifdef __OPENMP_BACKEND__    
            #pragma omp parallel
            {
                // get the thread private context
                // -> note that the context variable is thread private as it is defined inside the parallel section
                assert(static_cast<int>(threadPrivateContexts_.size()) > omp_get_thread_num());
                Context context = threadPrivateContexts_[omp_get_thread_num()];
#else
            {
                // get the thread private context 
                // -> note as we are compiling without OpenMP there is only one context object
                Context context = threadPrivateContexts_[0];
#endif

#ifndef NDEBUG
                // reset the bounding boxes of all iterators
                modify_tuple<
                    ResetBoundingBoxFunctor,
                    ParameterElementIndexes
                >(context.dataTuple());  
#endif

                // update the pointers
                modify_2_tuples<
                    UpdateApplyFunctor,
                    ParameterElementIndexes
                >(parameterTuple_, context.dataTuple());

                // init variables representing the memorized position of the context to the origin
                // note that these variables are thread private as they are defined inside the parallel section
                int currentI = 0;
                int currentJ = 0;
#ifdef __OPENMP_BACKEND__
                #pragma omp for nowait
#endif
                for(int i = 0;  i < numberOfBlocks; ++i)
                {              
                    // get the origin coordinates of the block
                    context.SetBlockConfiguration(launchConfiguration.blockConfigurations()[i]);
                    const int iOrigin = launchConfiguration.blockConfigurations()[i].iBlockIndex * BlockSize::ISize::value;
                    const int jOrigin = launchConfiguration.blockConfigurations()[i].jBlockIndex * BlockSize::JSize::value;
                    
                    // advance memorized positions in i if necessary
                    const int diffI = iOrigin - currentI;
                    currentI = iOrigin;
                    if(diffI != 0)
                    {
                        context.template AdvanceMemorizedPositionInI<DataFieldParameterIndexes>(diffI);
                    }
                    // advance memorized positions in j if necessary
                    const int diffJ = jOrigin - currentJ;
                    currentJ = jOrigin;
                    if(diffJ != 0)
                    {
                        context.template AdvanceMemorizedPositionInJ<DataFieldParameterIndexes>(diffJ);
                    }

                    // execute the stencil sweep groups
                    apply_to_all<
                        StencilSweepGroupFunctorOpenMP,
                        TStencilSweepGroupDescriptors,
                        Context
                    >(context);                        
                }
                
#ifndef NDEBUG
                // in debug mode finish all threads and compute the access bounding boxes
#ifdef __OPENMP_BACKEND__
                #pragma omp critical
#endif
                {
                    // merge the bounding boxes back
                    assert(threadPrivateContexts_.size() > 0);
                    modify_2_tuples<
                        MergeBoundingBoxFunctor,
                        DataFieldParameterIndexes
                    >(context.dataTuple(), threadPrivateContexts_[0].dataTuple());  
                }
#endif
            }
        }
    }

    // launch configuration manager
    LaunchConfigurationManager<LaunchConfiguration> launchConfigurationManager_;

    // contexts and input tuples
    std::vector<Context> threadPrivateContexts_; // keep context for every OpenMP thread
    ParameterTupleType parameterTuple_;
    std::vector<TemporaryFieldTupleType> threadPrivateTemporaryFieldTuples_;
};
  
// create an OpenMP stencil
template<
    typename TStencilConfiguration, 
    typename TParameterIJBoundary,
    typename TParameterTupleElements,
    typename TTemporaryFields,
    typename TStencilSweepGroupDescriptors>
struct create_stencil
{
    typedef StencilOpenMP<TStencilConfiguration, TParameterIJBoundary, TParameterTupleElements, TTemporaryFields, TStencilSweepGroupDescriptors> type;
};

