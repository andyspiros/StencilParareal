
// specialization terminating the modification recursion
template<bool VDone = true>
struct ModifyArrayImpl
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TData,
        typename TArray>
    MODIFY_ARRAY_QUALIFIER
    static void Execute(TArray& array, typename parameter_type<TData>::type data) {}

    // specialization without a data parameter
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TArray>
    MODIFY_ARRAY_QUALIFIER
    static void Execute(TArray& array) {}
};

// specialization recursively applying the functor to the array elements 
// addressed with TIndex
template<>
struct ModifyArrayImpl<false>
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TData,
        typename TArray>
    MODIFY_ARRAY_QUALIFIER 
    static void Execute(TArray& array, typename parameter_type<TData>::type data) 
    {
        // run the functor
        TFunctor::Do(
            array.At(static_cast<typename boost::mpl::deref<TIterator>::type*>(0)), 
            data
        );
        
        // recursively modify the other elements
        typedef typename boost::mpl::next<TIterator>::type Iter;
        ModifyArrayImpl< boost::is_same<Iter,TLastIterator>::value >::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor, 
            TData
        >(array, data);
    }

    // specialization without a data parameter
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TArray>
    MODIFY_ARRAY_QUALIFIER 
    static void Execute(TArray& array) 
    {
        // run the functor
        TFunctor::Do(
            array.At(static_cast<typename boost::mpl::deref<TIterator>::type*>(0))
        );

        // recursively modify the other elements
        typedef typename boost::mpl::next<TIterator>::type Iter;
        ModifyArrayImpl< boost::is_same<Iter,TLastIterator>::value >::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor
        >(array);
    }
};

/**
* Method modifying selected array elements
* @param array array to modify
* @param data data structure holding the state necessary for the modification
*/
template<   
    typename TFunctor,
    typename TIndexes,
    typename TData,
    typename TArray>
MODIFY_ARRAY_QUALIFIER
inline void modify_array(TArray& array, typename parameter_type<TData>::type data)
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TIndexes>::value);

    // note the array might be const
    BOOST_STATIC_ASSERT(is_array<typename boost::remove_const<TArray>::type>::value);

    typedef typename boost::mpl::begin<TIndexes>::type First;
    typedef typename boost::mpl::end<TIndexes>::type Last;

    // recursively apply functor to the data
    ModifyArrayImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor, 
        TData
    >(array, data);
}

/**
* Method modifying selected array elements (specialization without a data parameter)
* @param array array to modify
*/
template<   
    typename TFunctor,
    typename TIndexes,
    typename TArray>
MODIFY_ARRAY_QUALIFIER
inline void modify_array(TArray& array)
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TIndexes>::value);
    // note the array might be const
    BOOST_STATIC_ASSERT(is_array<typename boost::remove_const<TArray>::type>::value);

    typedef typename boost::mpl::begin<TIndexes>::type First;
    typedef typename boost::mpl::end<TIndexes>::type Last;

    // recursively apply functor to the data
    ModifyArrayImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor
    >(array);
}

