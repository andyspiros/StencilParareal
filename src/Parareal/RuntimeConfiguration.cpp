#include "RuntimeConfiguration.h"
#include "boost/program_options.hpp"

RuntimeConfiguration::RuntimeConfiguration(int argc, char **argv)
    : nu_(1.), cx_(1.), cy_(1.), cz_(1.)
    , gridSize_(32), endTime_(0.05)
{
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce this help message")
        ("nu", po::value<double>(), "Heat coefficient")
        ("cx", po::value<double>(), "Advection velocity in x direction")
        ("cy", po::value<double>(), "Advection velocity in y direction")
        ("cz", po::value<double>(), "Advection velocity in z direction")
        ("gridSize", po::value<int>(), "Number of grid points along each direction")
        ("endTime", po::value<double>(), "Time to simulate")
        ("cflFine", po::value<double>(), "CFL number of fine propagator")
        ("cflCoarse", po::value<double>(), "CFL number of coarse propagator")
        ("kmax", po::value<int>(), "Number of parareal iterations")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        std::exit(0);
    }

    if (vm.count("nu"))
    {
        this->set_nu(vm["nu"].as<double>());
    }

    if (vm.count("cx"))
    {
        this->set_cx(vm["cx"].as<double>());
    }

    if (vm.count("cy"))
    {
        this->set_cy(vm["cy"].as<double>());
    }

    if (vm.count("cz"))
    {
        this->set_cz(vm["cz"].as<double>());
    }

    if (vm.count("gridSize"))
    {
        this->set_gridSize(vm["gridSize"].as<int>());
    }

    if (vm.count("cflCoarse"))
    {
        this->set_cflCoarse(vm["cflCoarse"].as<double>());
    }

    if (vm.count("cflFine"))
    {
        this->set_cflFine(vm["cflFine"].as<double>());
    }

    if (vm.count("endTime"))
    {
        this->set_endTime(vm["endTime"].as<double>());
    }

    if (vm.count("kmax"))
    {
        this->set_endTime(vm["kmax"].as<int>());
    }

}
