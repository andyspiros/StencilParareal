
// tuple of size iteration index
template<typename TTupleElements>
class TupleImpl<TTupleElements, ITERATION_INDEX>
{
public:
    // define element types and indexes
    BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_DEFINE_INDEX_AND_ELEMENT_TYPES, void)

    TUPLE_QUALIFIER
    TupleImpl() {};
    TUPLE_QUALIFIER
    ~TupleImpl() {};

    TUPLE_QUALIFIER
    void Init(BOOST_PP_ENUM_BINARY_PARAMS(ITERATION_INDEX, T, val)) 
    {
        BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_INIT, void)
    };

    TUPLE_QUALIFIER
    TupleImpl(const TupleImpl& other) { *this = other; }
    TUPLE_QUALIFIER
    TupleImpl& operator= (const TupleImpl& other)
    {
        BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_ASSIGN, void)
        // by convention
        return *this;
    }

    // return void* if index is wrong
    template<typename T>
    TUPLE_QUALIFIER
    boost::mpl::void_ operator() (T*) { return boost::mpl::void_(); } 

    BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_OPERATOR, void)
    BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_CONST_OPERATOR, void)

private:
    BOOST_PP_REPEAT(ITERATION_INDEX, TEXT_MEMBER, void)
};


