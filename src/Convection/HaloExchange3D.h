#ifndef HALO_EXCHANGE_H_
#define HALO_EXCHANGE_H_

#include "GCL.h"
#include "utils/layout_map.h"
#include "utils/boollist.h"
#include "L2/include/halo_exchange.h"
#include "L2/include/descriptors.h"
#include "L3/include/proc_grids_3D.h"
#include "L3/include/Halo_Exchange_3D.h"

#include "SharedInfrastructure.h"

template<typename DataFieldType>
class HaloExchange3D
{
public:
    typedef typename DataFieldType::ValueType ValueType;

private:
    // calculate the index of StrideDim in the storage order
    template<typename StorageOrder, typename Dimension, Dimension VDimension>
    struct GetDimPosition
    {
        typedef typename boost::mpl::eval_if<
            boost::mpl::contains<
                StorageOrder,
                boost::mpl::integral_c<Dimension, VDimension>
            >,
            boost::mpl::integral_c<int,
                boost::mpl::distance<
                    typename boost::mpl::begin<StorageOrder>::type,
                    typename boost::mpl::find<
                        StorageOrder,
                        boost::mpl::integral_c<Dimension, VDimension>
                    >::type
                >::value
            >,
            boost::mpl::integral_c<int, -1>
        >::type Index;
    };

    typedef GCL::MPI_3D_process_grid_t<GCL::gcl_utils::boollist<3> > GridType;
    typedef typename DataFieldType::StorageFormat::StorageOrder StorageOrder;
    typedef typename GetDimPosition<StorageOrder, typeof(cDimI), cDimI>::Index IPos;
    typedef typename GetDimPosition<StorageOrder, typeof(cDimJ), cDimJ>::Index JPos;
    typedef typename GetDimPosition<StorageOrder, typeof(cDimK), cDimK>::Index KPos;
    typedef GCL::layout_map<IPos::value, JPos::value, KPos::value> FieldLayout;
    typedef typename GridType::period_type PeriodType;

#ifdef __CUDA_BACKEND__
    typedef GCL::gcl_gpu CPUorGPUType;
#else
    typedef GCL::gcl_cpu CPUorGPUType;
#endif

    typedef GCL::halo_exchange_dynamic_ut<
        FieldLayout,
        GCL::layout_map<0, 1, 2>,
        ValueType,
        3,
        CPUorGPUType,
        GCL::version_manual
    > HandlerType;

public:
    HaloExchange3D(bool periodicI, bool periodicJ, bool periodicK, MPI_Comm comm)
        : hd_(PeriodType(periodicI, periodicJ, periodicK), comm), setup_(false)
    {
    }

    void registerField(DataFieldType& field)
    {
        IJKBoundary boundary = field.boundary();
        IJKSize csize = field.calculationDomain();
        IJKSize psize = field.storage().paddedSize();
        IJKSize asize = field.storage().allocatedSize();

        hd_.template add_halo<0>(
                -boundary.iMinusOffset(),
                boundary.iPlusOffset(),
                -boundary.iMinusOffset(),
                -boundary.iMinusOffset() + csize.iSize() -1,
                psize.iSize()
            );

        hd_.template add_halo<1>(
                -boundary.jMinusOffset(),
                boundary.jPlusOffset(),
                -boundary.jMinusOffset(),
                -boundary.jMinusOffset() + csize.jSize() -1,
                psize.jSize()
            );

        hd_.template add_halo<2>(
                -boundary.kMinusOffset(),
                boundary.kPlusOffset(),
                -boundary.kMinusOffset(),
                -boundary.kMinusOffset() + csize.kSize() -1,
                psize.kSize()
            );

        fields_.push_back(&field);
    }

    void exchange()
    {
        // Setup
        if (!setup_)
        {
            hd_.setup(fields_.size());
            fieldStorages_.resize(fields_.size());
            setup_ = true;
        }

        // Fill vector with pointers
        for (int i = 0; i < fields_.size(); ++i)
        {
            IJKBoundary boundary = fields_[i]->boundary();
            fieldStorages_[i] = fields_[i]->storage().pStorageBase();
        }
        
        // Exchange
        hd_.pack(fieldStorages_);
        hd_.exchange();
        hd_.unpack(fieldStorages_);
    }
private:
    HandlerType hd_;
    bool setup_;
    std::vector<DataFieldType*> fields_;
    std::vector<ValueType*> fieldStorages_;
};


#endif // HALO_EXCHANGE_H

