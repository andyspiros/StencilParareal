#include <iostream>
#include <iomanip>

#include "GCL.h"
#include "utils/layout_map.h"
#include "utils/boollist.h"
#include "L2/include/halo_exchange.h"
#include "L2/include/descriptors.h"
#include "L3/include/proc_grids_3D.h"
#include "L3/include/Halo_Exchange_3D.h"

#include "SharedInfrastructure.h"

#include "MatFile.h"


template<typename DataFieldType, typename GridType>
class HaloExchange3D
{

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

    typedef typename DataFieldType::ValueType ValueType;
    typedef typename DataFieldType::StorageFormat::StorageOrder StorageOrder;
    typedef typename GetDimPosition<StorageOrder, typeof(cDimI), cDimI>::Index IPos;
    typedef typename GetDimPosition<StorageOrder, typeof(cDimJ), cDimJ>::Index JPos;
    typedef typename GetDimPosition<StorageOrder, typeof(cDimK), cDimK>::Index KPos;
    typedef GCL::layout_map<IPos::value, JPos::value, KPos::value> FieldLayout;
    typedef typename GridType::period_type PeriodType;

    typedef GCL::halo_exchange_dynamic_ut<
        FieldLayout,
        GCL::layout_map<0, 1, 2>,
        ValueType,
        3,
        GCL::gcl_cpu,
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

        ValueType* origin = &field(
                boundary.iMinusOffset(),
                boundary.jMinusOffset(),
                boundary.kMinusOffset()
            );
        fields_.push_back(origin);
    }

    void exchange()
    {
        if (!setup_)
            hd_.setup(fields_.size());
        hd_.pack(fields_);
        hd_.exchange();
        hd_.unpack(fields_);
    }
private:
    HandlerType hd_;
    bool setup_;
    std::vector<ValueType*> fields_;
};


template<typename DataFieldType>
void printField(const DataFieldType& field)
{
    IJKSize size = field.calculationDomain();
    IJKBoundary boundary = field.boundary();
    int iStart = boundary.iMinusOffset();
    int jStart = boundary.jMinusOffset();
    int kStart = boundary.kMinusOffset();
    int iEnd = size.iSize() + boundary.iPlusOffset();
    int jEnd = size.jSize() + boundary.jPlusOffset();
    int kEnd = size.kSize() + boundary.kPlusOffset();

    for (int k = kStart; k < kEnd; ++k)
    {
        for (int i = iStart; i < iEnd; ++i)
        {
            for (int j = jStart; j < jEnd; ++j)
            {
                std::cout << std::hex << (int)field(i, j, k) << "  ";
            }
            std::cout << "\n";
        }
        std::cout << "\n\n";
    }
}

int main()
{
    // typedef DataFieldOpenMP<
    //         Real,
    //         DataFieldStorageFormat<
    //             OpenMPIJBoundary,
    //             StorageOrder::JKI,
    //             DataFieldAlignment<cDimI, 1>
    //         >
    //     >
    //     MyRealField;
    typedef IJKRealField MyRealField;
    typedef GCL::gcl_utils::boollist<3> CyclicType;
    typedef GCL::MPI_3D_process_grid_t<CyclicType> GridType;

    const double pi = 3.14159265358979;

    // Initialize MPI
    GCL::GCL_Init();
    int isroot = GCL::PID == 0;

    int commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);

    int dims[3] = {0, 0, 1};
    int periods[3] = {1, 1, 1};
    MPI_Dims_create(GCL::PROCS, 3, dims);
    MPI_Comm comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 0, &comm);

    // Initialize GCL grid
    GridType grid(CyclicType(true, true, true), comm);
    int myPI, myPJ, myPK;
    grid.coords(myPI, myPJ, myPK);

    HaloExchange3D<MyRealField, GridType> he(true, true, true, comm);

    // Initializing field
    IJKSize domain; domain.Init(32, 32, 32);
    KBoundary kboundary; kboundary.Init(-3, 3);
    MyRealField field;
    field.Init("field", domain, kboundary);

    double dx = 1. / (domain.iSize()+1);

    for (int i = 0; i < domain.iSize(); ++i)
        for (int j = 0; j < domain.jSize(); ++j)
            for (int k = 0; k < domain.kSize(); ++k)
            {
                double x = (i+1) * dx;
                double y = (j+1) * dx;
                double z = (k+1) * dx;
                field(i, j, k) = sin(2.*pi*x) * sin(2.*pi*y) * sin(2.*pi*z);
            }

    he.registerField(field);

    MatFile mf("test.mat");

    mf.addField(field, 0);
    he.exchange();
    mf.addField(field, 1);

}

