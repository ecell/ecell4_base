template<typename Timpl_>
struct domain_id_pair_converter
{
    typedef Timpl_ native_type;

    struct to_python_converter
    {
        static PyObject* convert(native_type const& v)
        {
            return boost::python::incref(boost::python::make_tuple(
                v.first, v.second).ptr());
        }
    };

    static void __register()
    {
        boost::python::to_python_converter<native_type, to_python_converter>();
        peer::util::to_native_converter<native_type,
            peer::converters::pytuple_to_tuple_converter<native_type> >();
    }
};

class _EGFRDSimulator: public b::EGFRDSimulator, public boost::python::wrapper<b::EGFRDSimulator>
{
    typedef b::EGFRDSimulator base_type;
    typedef base_type::world_type world_type;
    typedef base_type::network_rules_type network_rules_type;
    typedef base_type::rng_type rng_type;
    typedef base_type::domain_id_pair_generator domain_id_pair_generator;

public:
    _EGFRDSimulator(world_type& world, rng_type& rng,
                    network_rules_type const& network_rules)
        : base_type(world, rng, network_rules) {}

    virtual void step()
    {
        get_override("step")();
    }
};

class_<_EGFRDSimulator, boost::noncopyable>("_EGFRDSimulator",
            init<CyclicWorld&,
                 GSLRandomNumberGenerator&,
                 NetworkRulesWrapper const&>())
        .def("new_spherical_shell", &_EGFRDSimulator::new_spherical_shell)
        .def("new_cylindrical_shell", &_EGFRDSimulator::new_spherical_shell)
        .def("get_spherical_shell",
            &_EGFRDSimulator::new_spherical_shell,
            return_value_policy<return_by_value>())
        .def("get_cylindrical_shell",
            &_EGFRDSimulator::new_spherical_shell,
            return_value_policy<return_by_value>())
        .def("get_domain", &_EGFRDSimulator::get_domain)
        .def("update_domain", &_EGFRDSimulator::update_domain)
        .def("__iter__", &_EGFRDSimulator::get_domains,
                return_value_policy<return_by_value>())
        ;

    peer::wrappers::generator_wrapper<ptr_generator<_EGFRDSimulator::domain_id_pair_generator, std::auto_ptr<_EGFRDSimulator::domain_id_pair_generator> > >::__register_class("DomainIDPairGenerator");

    domain_id_pair_converter<DomainID>::__register();

