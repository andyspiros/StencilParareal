
// specialization terminating the modification recursion
template<bool VDone = true>
struct ModifyTupleImpl
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TData,
        typename TTuple>
    MODIFY_TUPLE_QUALIFIER
    static void Execute(TTuple& tuple, typename parameter_type<TData>::type data) {}

    // specialization without a data parameter
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TTuple>
    MODIFY_TUPLE_QUALIFIER
    static void Execute(TTuple& tuple) {}
};

// specialization recursively applying the functor to the tuple elements 
// addressed with TIndex
template<>
struct ModifyTupleImpl<false>
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TData,
        typename TTuple>
    MODIFY_TUPLE_QUALIFIER 
    static void Execute(TTuple& tuple, typename parameter_type<TData>::type data) 
    {
        // run the functor
        TFunctor::Do(
            tuple(static_cast<typename boost::mpl::deref<TIterator>::type*>(0)), 
            data
        );
        
        // recursively modify the other elements
        typedef typename boost::mpl::next<TIterator>::type Iter;
        ModifyTupleImpl< boost::is_same<Iter,TLastIterator>::value >::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor, 
            TData
        >(tuple, data);
    }

    // specialization without a data parameter
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TTuple>
    MODIFY_TUPLE_QUALIFIER 
    static void Execute(TTuple& tuple) 
    {
        // run the functor
        TFunctor::Do(
            tuple(static_cast<typename boost::mpl::deref<TIterator>::type*>(0))
        );

        // recursively modify the other elements
        typedef typename boost::mpl::next<TIterator>::type Iter;
        ModifyTupleImpl< boost::is_same<Iter,TLastIterator>::value >::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor
        >(tuple);
    }
};

/**
* Method modifying selected tuple elements.
* @param tuple tuple to modify
* @param data data structure holding the state necessary for the modification
*/
template<   
    typename TFunctor,
    typename TIndexes,
    typename TData,
    typename TTuple>
MODIFY_TUPLE_QUALIFIER
inline void modify_tuple(TTuple& tuple, typename parameter_type<TData>::type data)
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TIndexes>::value);
    // note the tuple might be const
    BOOST_STATIC_ASSERT(is_tuple<typename boost::remove_const<TTuple>::type>::value);

    typedef typename boost::mpl::begin<TIndexes>::type First;
    typedef typename boost::mpl::end<TIndexes>::type Last;

    // recursively apply functor to the data
    ModifyTupleImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor, 
        TData
    >(tuple, data);
}

/**
* Method modifying selected tuple elements
* @param tuple tuple to modify
*/
template<   
    typename TFunctor,
    typename TIndexes,
    typename TTuple>
MODIFY_TUPLE_QUALIFIER
inline void modify_tuple(TTuple& tuple)
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TIndexes>::value);
    // note the tuple might be const
    BOOST_STATIC_ASSERT(is_tuple<typename boost::remove_const<TTuple>::type>::value);

    typedef typename boost::mpl::begin<TIndexes>::type First;
    typedef typename boost::mpl::end<TIndexes>::type Last;

    // recursively apply functor to the data
    ModifyTupleImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor
    >(tuple);
}

// specialization terminating the modification recursion
template<bool VDone = true>
struct Modify2TuplesImpl
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TData,
        typename TTuple1,
        typename TTuple2>
    MODIFY_TUPLE_QUALIFIER
    static void Execute(TTuple1& tuple1, TTuple2& tuple2, typename parameter_type<TData>::type data) {}

    // specialization without a data parameter
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TTuple1,
        typename TTuple2>
    MODIFY_TUPLE_QUALIFIER
    static void Execute(TTuple1& tuple1, TTuple2& tuple2) {}
};

// specialization recursively applying the functor to the tuple elements addressed with TIndex
template<>
struct Modify2TuplesImpl<false>
{
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TData,
        typename TTuple1,
        typename TTuple2>
    MODIFY_TUPLE_QUALIFIER
    static void Execute(TTuple1& tuple1, TTuple2& tuple2, typename parameter_type<TData>::type data) 
    {
        // run the functor
        TFunctor::Do(
            tuple1(static_cast<typename boost::mpl::deref<TIterator>::type*>(0)), 
            tuple2(static_cast<typename boost::mpl::deref<TIterator>::type*>(0)), 
            data
        );

        // recursively modify the other elements
        typedef typename boost::mpl::next<TIterator>::type Iter;
        Modify2TuplesImpl<boost::is_same<Iter,TLastIterator>::value>::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor, 
            TData
        >(tuple1, tuple2, data);
    }

    // specialization without a data parameter
    template<
        typename TIterator,
        typename TLastIterator,
        typename TFunctor,
        typename TTuple1,
        typename TTuple2>
    MODIFY_TUPLE_QUALIFIER
    static void Execute(TTuple1& tuple1, TTuple2& tuple2) 
    {
        // run the functor
        TFunctor::Do(
            tuple1(static_cast<typename boost::mpl::deref<TIterator>::type*>(0)), 
            tuple2(static_cast<typename boost::mpl::deref<TIterator>::type*>(0))
        );

        // recursively modify the other elements
        typedef typename boost::mpl::next<TIterator>::type Iter;
        Modify2TuplesImpl<boost::is_same<Iter,TLastIterator>::value>::template Execute<
            Iter, 
            TLastIterator, 
            TFunctor
        >(tuple1, tuple2);
    }
};

/**
* Method modifying selected tuple elements of two tuples. Normally used to convert one tuple to another one
* @param tuple1 first tuple to modify
* @param tuple1 second tuple to modify
* @param data data structure holding the state necessary for the modification
*/
template<   
    typename TFunctor,
    typename TIndexes,
    typename TData,
    typename TTuple1,
    typename TTuple2>
MODIFY_TUPLE_QUALIFIER
inline void modify_2_tuples(TTuple1& tuple1, TTuple2& tuple2, typename parameter_type<TData>::type data)
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TIndexes>::value);
    // note the tuples might be const
    BOOST_STATIC_ASSERT(is_tuple<typename boost::remove_const<TTuple1>::type>::value);
    BOOST_STATIC_ASSERT(is_tuple<typename boost::remove_const<TTuple2>::type>::value);

    typedef typename boost::mpl::begin<TIndexes>::type First;
    typedef typename boost::mpl::end<TIndexes>::type Last;

    // recursively apply functor to the data
    Modify2TuplesImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor, 
        TData
    >(tuple1, tuple2, data);
}

/**
* Method modifying selected tuple elements of two tuples
* (note that the method is usually used in order to covert a tuple)
* @param tuple1 first tuple to modify
* @param tuple1 second tuple to modify
*/
template<   
    typename TFunctor,
    typename TIndexes,
    typename TTuple1,
    typename TTuple2>
MODIFY_TUPLE_QUALIFIER
inline void modify_2_tuples(TTuple1& tuple1, TTuple2& tuple2)
{
    BOOST_STATIC_ASSERT(boost::mpl::is_sequence<TIndexes>::value);
    // note the tuples might be const
    BOOST_STATIC_ASSERT(is_tuple<typename boost::remove_const<TTuple1>::type>::value);
    BOOST_STATIC_ASSERT(is_tuple<typename boost::remove_const<TTuple2>::type>::value);

    typedef typename boost::mpl::begin<TIndexes>::type First;
    typedef typename boost::mpl::end<TIndexes>::type Last;

    // recursively apply functor to the data
    Modify2TuplesImpl<boost::is_same<First,Last>::value>::template Execute<
        First, 
        Last, 
        TFunctor
    >(tuple1, tuple2);
}

  
