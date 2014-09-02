#include <iostream>
#include <fstream>
#include <iomanip>
#include "RuntimeConfiguration.h"
#include "ConvectionSolution.h"
#include "mpi.h"
#include "Convection.h"
#include "SharedInfrastructure.h"
#include "StencilFramework.h"
#include "MatFile.h"
#include "Parareal.h"

// Power constants (Watt)
const double powerNetwork = 25.;
const double powerBlower = 14.21;

/**
 * Synchronzes the host storage of the field in the case of a GPU build
 */
template<typename TDataField>
inline void SynchronizeHost(TDataField& field)
{
#ifdef __CUDA_BACKEND__
    field.SynchronizeHostStorage();
#endif
}

/**
 * Synchronzes the device storage of the field in the case of a GPU build
 */
template<typename TDataField>
inline void SynchronizeDevice(TDataField& field)
{
#ifdef __CUDA_BACKEND__
    field.SynchronizeDeviceStorage();
#endif
}

/**
 * Synchronizes the CUDA device in the case of a GPU build
 */
inline void SynchronizeCUDA()
{
#ifdef __CUDA_BACKEND__
    cudaDeviceSynchronize();
#endif
}

/**
 * Read power management file and return contained value
 */
double read_pm_file(const std::string &fname) {
    double result = 0.;
    std::ifstream fid(fname.c_str());
    fid >> result;
    return result;
}

/**
 * Get node energy
 */
double energy()
{
    return read_pm_file("/sys/cray/pm_counters/energy");
}

/**
 * Get device energy
 */
double deviceEnergy()
{
    return read_pm_file("/sys/cray/pm_counters/accel_energy");
}

/**
 * Get total energy, single-process case
 */
double totalEnergy(double energyStart, double energyEnd)
{
    return energyEnd - energyStart;
}

/**
 * Get total energy, distributed case
 */
double totalEnergy(double energyStart, double energyEnd, MPI_Comm comm)
{
    double myEnergy = energyEnd - energyStart;
    double totEnergy;
    MPI_Allreduce(&myEnergy, &totEnergy, 1, MPI_DOUBLE, MPI_SUM, comm);
    return totEnergy;
}

/**
 * Gives the infinity norm of the error of the provided field
 * by checking against the reference solution
 */
template<typename TDataField>
double computeErrorReference(const TDataField& field, const TDataField& reference)
{
#ifdef __CUDA_BACKEND__
    const bool synchronizeField = field.isDeviceUpToDate();
    const bool synchronizeReference = reference.isDeviceUpToDate();
    if (synchronizeField)
        field.SynchronizeHostStorage();
    if (synchronizeReference)
        reference.SynchronizeHostStorage();
#endif

    const IJKSize& size = field.calculationDomain();
    const int isize = size.iSize();
    const int jsize = size.jSize();
    const int ksize = size.kSize();

    const IJKBoundary boundary = field.boundary();
    const int istart = -boundary.iMinusOffset();
    const int jstart = -boundary.jMinusOffset();
    const int kstart = -boundary.kMinusOffset();

    double error = 0., maxfield = 0.;

    for (int i = 0; i < isize; ++i)
        for (int j = 0; j < jsize; ++j)
            for (int k = 0; k < ksize; ++k)
            {
                double e = std::abs(field(i,j,k) - reference(i,j,k));
                error = std::max(error, e);
                maxfield = std::max(std::abs(field(i, j, k)), maxfield);
            }

#ifdef __CUDA_BACKEND__
    if (synchronizeField)
        field.SynchronizeDeviceStorage();
    if (synchronizeReference)
        reference.SynchronizeDeviceStorage();
#endif

    return error/maxfield;
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int commsize, commrank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
    const bool isRoot = commrank == 0;
    const bool isLast = commrank == commsize - 1;

    if (isRoot)
        std::cout << "Initialization...\n" << std::endl;

    RuntimeConfiguration conf(argc, argv);

    // Compute my start and end time
    const double timeStart = conf.timeSliceSize() * commrank;
    const double timeEnd = conf.timeSliceSize() * (commrank + 1);

    if (isRoot)
    std::cout << "Running with:\n"
        << " - initial diffusion coefficient: " << conf.nu0() << "\n"
        << " - frequence of diffusion coefficient: " << conf.nufreq() << "\n"
        << " - advection velocity in x: " << conf.cx() << "\n"
        << " - advection velocity in y: " << conf.cy() << "\n"
        << " - advection velocity in z: " << conf.cz() << "\n"
        << " - spatial discretization step: " << conf.dx() << "\n"
        << " - endtime: " << conf.endTime() << "\n"
        << " - number of time slices: " << conf.timeSlices() << "\n"
        << " - time slice size: " << conf.timeSliceSize() << "\n"
        << " - CFL fine: " << conf.cflFine() << "\n"
        << " - CFL coarse: " << conf.cflCoarse() << "\n"
        << " - timestep size fine: " << conf.dtFine() << "\n"
        << " - timestep size coarse: " << conf.dtCoarse() << "\n"
        << " - timesteps per slice fine propagator: " << conf.timeStepsFinePerTimeSlice() << "\n"
        << " - timesteps per slice coarse propagator: " << conf.timeStepsCoarsePerTimeSlice() << "\n"
        << " - parareal iterations: " << conf.kmax() << "\n"
        << " - asynchronous communications: " << (conf.async() ? "Enabled" : "Disabled") << "\n"
        << " - intermediate fields in mat files: " << (conf.mat() ? "Yes" : "No") << "\n"
        << std::endl;

    // Calculation domain and boundaries
    IJKSize domain; domain.Init(conf.gridSize(), conf.gridSize(), conf.gridSize());
    KBoundary kboundary; kboundary.Init(-convectionBoundaryLines, convectionBoundaryLines);

    // Initialize fields
    ConvectionField q, qinitial;
    q.Init("q", domain, kboundary);
    qinitial.Init("qinitial", domain, kboundary);
    Convection convection(conf.gridSize(), conf.gridSize(), conf.gridSize(), conf.dx(), conf.nu0(), conf.nufreq(), conf.cx(), conf.cy(), conf.cz());

    // Initialize parareal
    Parareal<Convection, ConvectionField> parareal(convection, qinitial, q, timeStart, conf, MPI_COMM_WORLD);

    if (conf.mode() == ModeCompare)
    {
        // Measure time required by convection
        const int tauSamples = 4;
        double tauF = MPI_Wtime();
        convection.DoRK4(qinitial, qinitial, 0., conf.dtFine(), tauSamples*conf.timeStepsFinePerTimeSlice());
        SynchronizeCUDA();
        tauF = MPI_Wtime() - tauF;
        double tauG = MPI_Wtime();
        convection.DoEuler(qinitial, qinitial, 0., conf.dtCoarse(), tauSamples*conf.timeStepsCoarsePerTimeSlice());
        SynchronizeCUDA();
        tauG = MPI_Wtime() - tauG;

        const double tauRatio = tauG / tauF;
        const double Nit_Np = static_cast<double>(conf.kmax()) / commsize;
        const double maxSpeedup = 1. / (tauRatio * (1. + Nit_Np) + Nit_Np);

        // Fill initial solution
        SynchronizeHost(qinitial);
        fillQ(qinitial, conf.nu0(), conf.nufreq(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
        SynchronizeDevice(qinitial);

        // Run serial
        MPI_Barrier(MPI_COMM_WORLD);
        double eserial = MPI_Wtime();
        parareal.DoSerial();
        eserial = MPI_Wtime() - eserial;

        // Save reference
        ConvectionField qreference = q;
        SynchronizeHost(qreference);

        // Fill initial solution
        SynchronizeHost(qinitial);
        fillQ(qinitial, conf.nu0(), conf.nufreq(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
        SynchronizeDevice(qinitial);

        // Run serial
        MPI_Barrier(MPI_COMM_WORLD);
        double eparallel = MPI_Wtime();
        parareal.DoParallel();
        eparallel = MPI_Wtime() - eparallel;

        // Output
        MPI_Barrier(MPI_COMM_WORLD);
        if (isLast)
        {
            double e = computeErrorReference(q, qreference);
            std::cout << "\n"
                << "Serial run time: " << eserial << "\n"
                << "Parallel run time: " << eparallel << "\n"
                << "Speedup: " << eserial / eparallel << "\n"
                << "Maximal speedup: " << maxSpeedup << "\n"
                << "Error at end: " << e << "\n"
                << std::endl;

            MatFile matfile("result.mat");
            matfile.addField("q", q);
            matfile.addField("qreference", qreference);
        }
    }
    else if (conf.mode() == ModeSerial)
    {
        // Fill initial solution
        SynchronizeHost(qinitial);
        fillQ(qinitial, conf.nu0(), conf.nufreq(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
        SynchronizeDevice(qinitial);

        // Run serial
        double e = MPI_Wtime();
        double energyStart = energy();
        double deviceEnergyStart = deviceEnergy();

        parareal.DoSerial();

        e = MPI_Wtime() - e;
        double energyEnd = energy();
        double deviceEnergyEnd = deviceEnergy();

        const double totDevice = totalEnergy(deviceEnergyStart, deviceEnergyEnd);
        const double totNode = totalEnergy(energyStart, energyEnd) - totDevice;
        const double totNetwork = e * powerNetwork;
        const double totBlower = e * powerBlower;
        const double totEnergy = totNode + totNetwork + totBlower;

        // Output
        MPI_Barrier(MPI_COMM_WORLD);
        if (isLast)
        {
            std::cout << "\n" << "Serial run time: " << e << "\n";
            std::printf("Node energy   : %8f J  (%8.3e W/node)\n", totNode   , totNode/e);
            std::printf("Device energy : %8f J  (%8.3e W/node)\n", totDevice , totDevice/e);
            std::printf("Network energy: %8f J  (%8.3e W/node)\n", totNetwork, totNetwork/e);
            std::printf("Blower energy : %8f J  (%8.3e W/node)\n", totBlower , totBlower/e);
            std::printf("Total energy  : %8f J  (%8.3e W/node)\n", totEnergy , totEnergy/e);
            std::cout << std::endl;
        }
    }
    else if (conf.mode() == ModeParallel)
    {
        // Fill initial solution
        SynchronizeHost(qinitial);
        fillQ(qinitial, conf.nu0(), conf.nufreq(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
        SynchronizeDevice(qinitial);

        // Run parallel
        MPI_Barrier(MPI_COMM_WORLD);
        double e = MPI_Wtime();
        double energyStart = energy();
        double deviceEnergyStart = deviceEnergy();

        parareal.DoParallel();

        MPI_Barrier(MPI_COMM_WORLD);
        e = MPI_Wtime() - e;
        double energyEnd = energy();
        double deviceEnergyEnd = deviceEnergy();

        const double totDevice = totalEnergy(deviceEnergyStart, deviceEnergyEnd, MPI_COMM_WORLD);
        const double totNode = totalEnergy(energyStart, energyEnd, MPI_COMM_WORLD) - totDevice;
        const double totNetwork = e * powerNetwork * commsize;
        const double totBlower = e * powerBlower * commsize;
        const double totEnergy = totNode + totNetwork + totBlower;

        // Output
        MPI_Barrier(MPI_COMM_WORLD);
        if (isLast)
        {
            std::cout << "\n" << "Parallel run time: " << e << "\n";
            std::printf("Node energy   : %8f J  (%8.3e W/node)\n", totNode   , totNode/e/commsize);
            std::printf("Device energy : %8f J  (%8.3e W/node)\n", totDevice , totDevice/e/commsize);
            std::printf("Network energy: %8f J  (%8.3e W/node)\n", totNetwork, totNetwork/e/commsize);
            std::printf("Blower energy : %8f J  (%8.3e W/node)\n", totBlower , totBlower/e/commsize);
            std::printf("Total energy  : %8f J  (%8.3e W/node)\n", totEnergy , totEnergy/e/commsize);
            std::cout << std::endl;
        }
    }
    else if (conf.mode() == ModeTiming)
    {
        // Fill initial solution
        SynchronizeHost(qinitial);
        fillQ(qinitial, conf.nu0(), conf.nufreq(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
        SynchronizeDevice(qinitial);

        // Run serial
        std::vector<double> times;
        MPI_Barrier(MPI_COMM_WORLD);
        parareal.DoTimedParallel(times);

        // Gather on root
        const int s = times.size();
        std::vector<double> timesGlobal;
        timesGlobal.resize(s * commsize);
        MPI_Gather(&times[0], s, MPI_DOUBLE, &timesGlobal[0], s, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Output
        if (isRoot)
        {
            std::cout << "\nTimes:\n";
            for(int i = 0; i < s; ++i)
            {
                for(int p = 0; p < commsize; ++p)
                {
                    std::cout << std::scientific << std::setprecision(6) << timesGlobal[p*s + i] << "   ";
                }
                std::cout << "\n";
            }
        }
    }

    // Finalize
    MPI_Finalize();

    return 0;
}

