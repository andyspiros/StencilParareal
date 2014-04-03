
// specialization terminating the recursion
template<bool VDone = true>
struct ApplyToAllImpl
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TApplyData1>
    APPLY_TO_ALL_QUALIFIER
    static void Execute(typename parameter_type<TApplyData1>::type) {}

    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TApplyData1,
        typename TApplyData2>
    APPLY_TO_ALL_QUALIFIER
    static void Execute(typename parameter_type<TApplyData1>::type, typename parameter_type<TApplyData2>::type) {}
};

// specialization recursively running the functor passing a data parameter
template<>
struct ApplyToAllImpl<false>
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TApplyData1>
    APPLY_TO_ALL_QUALIFIER
    static void Execute(typename parameter_type<TApplyData1>::type data1) 
    {
        // run the functor
        TFunctor::template Do<
            TApplyData1, typename boost::mpl::deref<TIterator>::type
        >(data1);

        // recursively iterator over the type list
        typedef typename boost::mpl::next<TIterator>::type Iter;
        ApplyToAllImpl<boost::is_same<Iter,TLastIterator>::value>::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor, 
            TApplyData1
        >(data1);
    }

    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TApplyData1,
        typename TApplyData2>
    APPLY_TO_ALL_QUALIFIER
    static void Execute(typename parameter_type<TApplyData1>::type data1, typename parameter_type<TApplyData2>::type data2) 
    {
        // run the functor
        TFunctor::template Do<
            TApplyData1, TApplyData2, typename boost::mpl::deref<TIterator>::type
        >(data1, data2);

        // recursively iterator over the type list
        typedef typename boost::mpl::next<TIterator>::type Iter;
        ApplyToAllImpl<boost::is_same<Iter,TLastIterator>::value>::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor, 
            TApplyData1,
            TApplyData2
        >(data1, data2);
    }
};

/**
* Method iterating over a type list and applying a functor to the elements of the type list
* @param data1 functor parameter
*/
template<
    typename TFunctor,
    typename TVector,
    typename TApplyData1>
APPLY_TO_ALL_QUALIFIER
inline void apply_to_all(typename parameter_type<TApplyData1>::type data1) 
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TVector>::value);

    typedef typename boost::mpl::begin<TVector>::type First;
    typedef typename boost::mpl::end<TVector>::type Last;

    ApplyToAllImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor, 
        TApplyData1
    >(data1);
}

/**
* Method iterating over a type list and applying a functor to the elements of the type list
* @param data1 functor parameter
* @param data2 functor parameter
*/
template<
    typename TFunctor,
    typename TVector,
    typename TApplyData1,
    typename TApplyData2>
APPLY_TO_ALL_QUALIFIER
inline void apply_to_all(typename parameter_type<TApplyData1>::type data1, typename parameter_type<TApplyData2>::type data2) 
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TVector>::value);

    typedef typename boost::mpl::begin<TVector>::type First;
    typedef typename boost::mpl::end<TVector>::type Last;

    ApplyToAllImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor, 
        TApplyData1,
        TApplyData2
    >(data1, data2);
}  
