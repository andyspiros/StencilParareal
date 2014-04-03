#ifndef BOOST_PP_IS_ITERATING

    #ifndef DATA_FIELD_STORAGE_STRIDES_INCLUDED
    #define DATA_FIELD_STORAGE_STRIDES_INCLUDED

    #include <boost/preprocessor/repetition.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/arithmetic/sub.hpp>
    #include <boost/preprocessor/arithmetic/add.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include <boost/mpl/integral_c.hpp>
    #include <boost/mpl/size.hpp>
    #include <boost/mpl/at.hpp>
    #include <boost/mpl/begin.hpp>
    #include <boost/mpl/contains.hpp>
    #include <boost/mpl/distance.hpp>
    #include <boost/mpl/find.hpp>
    #include <boost/mpl/eval_if.hpp>
    #include "Definitions.h"
    #include "Enums.h"
    #include "IJKSize.h"
    #include "DataFieldStorageFormat.h"

    // define the implementation class
    template<typename TStorageOrder, int VRank>
    class DataFieldStorageStridesImpl;

    /**
    * @class DataFieldStrides
    * Class handling the stride computation for data fields of various rank and storage order.
    */
    template<typename TStorageOrder>
    class DataFieldStorageStrides : public DataFieldStorageStridesImpl<TStorageOrder, boost::mpl::size<TStorageOrder>::value> {};

    // generate specializations
    #define BOOST_PP_ITERATION_PARAMS_1 (3, (0, MAX_FIELD_DIMENSION, "DataFieldStorageStrides.h"))
    #include BOOST_PP_ITERATE()

    #endif // DATA_FIELD_STORAGE_STRIDES_INCLUDED

#else // BOOST_PP_IS_ITERATING
    #define ITERATION_INDEX BOOST_PP_ITERATION()

    #define TEXT_CONSTRUCTOR(z, n, data) \
        BOOST_PP_CAT(BOOST_PP_CAT(stride, BOOST_PP_ADD(n, 1)), _) = -1;

    #define TEXT_ASSIGN(z, n, data) \
        BOOST_PP_CAT(BOOST_PP_CAT(stride, BOOST_PP_ADD(n, 1)), _) = \
        BOOST_PP_CAT(BOOST_PP_CAT(other.stride, BOOST_PP_ADD(n, 1)), _);

    #define TEXT_INIT(z, n, data) \
        BOOST_PP_CAT(BOOST_PP_CAT(stride, BOOST_PP_ADD(n, 1)), _) = \
        BOOST_PP_CAT(BOOST_PP_CAT(stride, n), _) * \
        sizes[boost::mpl::deref< BOOST_PP_CAT(Iter, n) >::type::value]; \
        typedef typename boost::mpl::next< BOOST_PP_CAT(Iter, n) >::type BOOST_PP_CAT(Iter, BOOST_PP_ADD(n, 1));

    #define TEXT_STRIDE(z, n, data) \
        __ACC_CPU__ \
        int strideByStorageOrderImpl(boost::mpl::integral_c<int, BOOST_PP_ADD(n, 1)>*) const \
        { return BOOST_PP_CAT(BOOST_PP_CAT(stride, BOOST_PP_ADD(n, 1)), _); }

    #define TEXT_MEMBER(z, n, data) \
        BOOST_PP_CAT(BOOST_PP_CAT(int stride, BOOST_PP_ADD(n, 1)), _);

    // define strides specialization
    template<typename TStorageOrder>
    class DataFieldStorageStridesImpl<TStorageOrder, ITERATION_INDEX>
    {
    public:
        typedef TStorageOrder StorageOrder;

        __ACC_CPU__
        DataFieldStorageStridesImpl() 
        {
            BOOST_PP_REPEAT(BOOST_PP_SUB(ITERATION_INDEX, 1), TEXT_CONSTRUCTOR, void)
        }
        __ACC_CPU__
        ~DataFieldStorageStridesImpl() {}

        __ACC_CPU__
        DataFieldStorageStridesImpl(const DataFieldStorageStridesImpl& other) { *this = other; }
        __ACC_CPU__
        DataFieldStorageStridesImpl& operator= (const DataFieldStorageStridesImpl& other) 
        {
            BOOST_PP_REPEAT(BOOST_PP_SUB(ITERATION_INDEX, 1), TEXT_ASSIGN, void)
            // by convention
            return *this;
        } 

        /**
        * Init sets up the strides.
        * @param size size in all dimensions
        */
        __ACC_CPU__
        void Init(const IJKSize& size)
        {
#if ITERATION_INDEX > 1
            // init the size structures
            int sizes[3];
            sizes[0] = size.iSize();
            sizes[1] = size.jSize();
            sizes[2] = size.kSize();
            
            // start iteration through storage orders
            int stride0_ = 1;
            typedef typename boost::mpl::begin<StorageOrder>::type Iter0;
            BOOST_PP_REPEAT(BOOST_PP_SUB(ITERATION_INDEX, 1), TEXT_INIT, void)
#endif
        }

        /**
        * ComputeStride the stride for a given index
        * @param i index in i direction
        * @param j index in j direction
        * @param k index in k direction
        */
        __ACC_CPU__
        int ComputeStride(const int i, const int j, const int k) const
        {
            return 
                i * strideByDimension<cDimI>() + 
                j * strideByDimension<cDimJ>() + 
                k * strideByDimension<cDimK>();
        }

    private:
        // different implementations for stride calculation
        template<Dimension VDimension>
        __ACC_CPU__
        int strideByDimension() const
        {
            // calculate the index of StrideDim in the storage order
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

            return strideByStorageOrderImpl(static_cast<Index*>(0));
        }

        // by default return stride 0
        template<int VIndex>
        __ACC_CPU__
        int strideByStorageOrderImpl(boost::mpl::integral_c<int, VIndex>*) const { return 0; }

#if ITERATION_INDEX > 0
        // define const stride and storage stride access overloads
        __ACC_CPU__ 
        int strideByStorageOrderImpl(boost::mpl::integral_c<int, 0>*) const { return 1; }
#endif

        BOOST_PP_REPEAT(BOOST_PP_SUB(ITERATION_INDEX, 1), TEXT_STRIDE, void)
        BOOST_PP_REPEAT(BOOST_PP_SUB(ITERATION_INDEX, 1), TEXT_MEMBER, void)
    };

    #undef ITERATION_INDEX
    #undef TEXT_CONSTRUCTOR
    #undef TEXT_ASSIGN
    #undef TEXT_INIT
    #undef TEXT_STRIDE
    #undef TEXT_MEMBER

#endif // BOOST_PP_IS_ITERATING


