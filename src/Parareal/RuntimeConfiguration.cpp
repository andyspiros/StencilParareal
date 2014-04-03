#include "RuntimeConfiguration.h"
#include "boost/program_options.hpp"
#include "mpi.h"
#include <string>
#include <algorithm>

RuntimeConfiguration::RuntimeConfiguration(int argc, char **argv)
    : nu0_(1.), nufreq_(0.), cx_(1.), cy_(1.), cz_(1.)
    , gridSize_(32), endTime_(0.05)
    , timeStepsFine_(128), timeStepsCoarse_(32)
    , kmax_(5), mat_(false), async_(true)
    , mode_(ModeCompare)
{
    // Retrieve MPI information
    MPI_Comm_size(MPI_COMM_WORLD, &timeSlices_);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const bool isRoot = !rank;

    // Create description of allowed options
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "Produce this help message")
        ("nu0", po::value<double>(), "Initial value for diffusion coefficient")
        ("nufreq", po::value<double>(), "Frequency for diffusion coefficient")
        ("cx", po::value<double>(), "Advection velocity in x direction")
        ("cy", po::value<double>(), "Advection velocity in y direction")
        ("cz", po::value<double>(), "Advection velocity in z direction")
        ("gridSize", po::value<int>(), "Number of grid points along each direction")
        ("endTime", po::value<double>(), "Time to simulate")
        ("timeStepsFine", po::value<int>(), "Total number of fine timesteps")
        ("timeStepsCoarse", po::value<int>(), "Total number of coarse timesteps")
        ("kmax", po::value<int>(), "Number of parareal iterations")
        ("mat", "Serialize intermediate fields")
        ("noasync", "Avoid asynchronous communication")
        ("mode", po::value<std::string>(), "Select among Compare, Serial and Parallel")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Parse options

    if (vm.count("help"))
    {
        std::cout << desc << "\n";
        std::exit(0);
    }

    if (vm.count("nu0"))
    {
        this->set_nu0(vm["nu0"].as<double>());
    }

    if (vm.count("nufreq"))
    {
        this->set_nufreq(vm["nufreq"].as<double>());
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

    if (vm.count("endTime"))
    {
        this->set_endTime(vm["endTime"].as<double>());
    }

    if (vm.count("timeStepsCoarse"))
    {
        this->set_timeStepsCoarse(vm["timeStepsCoarse"].as<int>());
    }

    if (vm.count("timeStepsFine"))
    {
        this->set_timeStepsFine(vm["timeStepsFine"].as<int>());
    }

    if (vm.count("kmax"))
    {
        this->set_kmax(vm["kmax"].as<int>());
    }

    if (vm.count("mat"))
    {
        this->set_mat(true);
    }

    if (vm.count("noasync"))
    {
        this->set_async(false);
    }

    if (vm.count("mode"))
    {
        std::string opt = vm["mode"].as<std::string>();
        std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
        if (opt == "compare")
            this->set_mode(ModeCompare);
        else if (opt == "parallel")
            this->set_mode(ModeParallel);
        else if (opt == "serial")
            this->set_mode(ModeSerial);
        else if (opt == "timing")
            this->set_mode(ModeTiming);
        else
        {
            std::cerr << "Mode " << vm["mode"].as<std::string>() << " not recognized.\n"
                << "Defaulting to Compare\n";
            this->set_mode(ModeCompare);
        }
    }

    if (timeStepsFine_ % timeSlices_ != 0)
    {
        timeStepsFine_ = (timeStepsFine_ / timeSlices_ + 1) * timeSlices_;
        if (isRoot)
        {
            std::cout << "Warning: timeStepsFine is not a multiple of timeSlices\n"
                << " -- Changing timeStepsFine to " << timeStepsFine_ << "\n";
        }
    }

    if (timeStepsCoarse_ % timeSlices_ != 0)
    {
        timeStepsCoarse_ = (timeStepsCoarse_ / timeSlices_ + 1) * timeSlices_;
        if (isRoot)
        {
            std::cout << "Warning: timeStepsCoarse is not a multiple of timeSlices\n"
                << " -- Changing timeStepsCoarse to " << timeStepsCoarse_ << "\n\n";
        }
    }
}

