
// array of size iteration index
template<typename TValue>
class Array<TValue, ITERATION_INDEX>
{
public:
    typedef typename boost::remove_const<TValue>::type ValueType;

#if ITERATION_INDEX > 1
    typedef BOOST_PP_CAT(BOOST_PP_CAT(boost::mpl::vector, ITERATION_INDEX), _c)<int, BOOST_PP_ENUM(ITERATION_INDEX, TEXT_NUMBER, void)> ElementIndexes;
#else
    typedef boost::mpl::vector0_c<int> ElementIndexes;
#endif
  
    ARRAY_QUALIFIER
    Array() {};
    ARRAY_QUALIFIER
    ~Array() {};

    ARRAY_QUALIFIER
    Array(const Array& other) { *this = other; }
    ARRAY_QUALIFIER
    Array& operator= (const Array& other)
    {
        BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_ASSIGN, void)
        // by convention
        return *this;
    }

    // return void for unknown indexes 
    template<int VIndex>
    ARRAY_QUALIFIER
    boost::mpl::void_ At(boost::mpl::integral_c<int, VIndex>*) { return boost::mpl::void_(); }

    BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_AT, void)
    BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_CONST_AT, void)

private:
    BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_MEMBER, void)
};


