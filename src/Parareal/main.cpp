#include <iostream>
#include <iomanip>
#include "RuntimeConfiguration.h"
#include "ConvectionSolution.h"
#include "mpi.h"
#include "Convection.h"
#include "SharedInfrastructure.h"
#include "StencilFramework.h"
#include "MatFile.h"

// Compile-time checks
#ifdef ENABLE_INTERNAL_TIMING
# ifdef ENABLE_ERROR
#  warning Enabling both internal timing and error computation after every iteration. Timings will be inaccurate
# endif
#endif

enum
{
    q_, qfine_, qcoarsenew_, qcoarseold_
};

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
 * Writes information about the internal timing if requested at compile-time
 */
inline void printTimer(double e, const char* str, std::ofstream& fs)
{
#ifdef ENABLE_INTERNAL_TIMING
    fs << " - " << str << " done in " << e << " msec\n";
#endif
}

template<typename TEnv>
struct UpdateStage
{
    STENCIL_STAGE(TEnv)

    STAGE_PARAMETER(FullDomain, q_)
    STAGE_PARAMETER(FullDomain, qfine_)
    STAGE_PARAMETER(FullDomain, qcoarseold_)
    STAGE_PARAMETER(FullDomain, qcoarsenew_)

    __ACC__
    static void Do(Context ctx, FullDomain)
    {
        ctx[q_::Center()] = 
            ctx[qcoarsenew_::Center()] + ctx[qfine_::Center()] - ctx[qcoarseold_::Center()];
    }
};

class ExternalTimer
{
public:
    void start()
    {
        SynchronizeCUDA();
        e_ = MPI_Wtime();
    }

    double elapsed()
    {
        SynchronizeCUDA();
        return (MPI_Wtime() - e_) * 1000.;
    }

private:
        double e_;
};

/**
 * Enabled only in case ENABLE_INTERNAL_TIMING is set at compile-time
 */
class InternalTimer
{
public:
    void start()
    {
#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        e_ = MPI_Wtime();
#endif
    }

    double elapsed()
    {
#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        return (MPI_Wtime() - e_) * 1000.;
#else
        return 0.;
#endif
    }

private:
        double e_;

};

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
    const int iend = istart + isize;
    const int jend = jstart + jsize;
    const int kend = kstart + ksize;

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
        << " - diffusion coefficient: " << conf.nu() << "\n"
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

    // Create log file
    std::ostringstream logfname;
    logfname << "Parareal_" << commrank << ".log";
    std::ofstream logfile(logfname.str().c_str());

    logfile << "Rank " << commrank << " integrates from " << timeStart
        << " to " << timeEnd << "\n" << std::endl;

    // Create mat file
    logfname.seekp(0);
    logfname << "Parareal_" << commrank << ".mat";
    MatFile matfile(logfname.str());

    // Initialize fields
    ConvectionField q, qinitial, qcoarseold, qcoarsenew, qfine, qreference;
    q.Init("q", domain, kboundary);
    qinitial.Init("qinitial", domain, kboundary);
    qcoarseold.Init("qcoarseold", domain, kboundary);
    qcoarsenew.Init("qcoarsenew", domain, kboundary);
    qfine.Init("qfine", domain, kboundary);
    qreference.Init("qreference", domain, kboundary);

    // Initialize update stencil
    Stencil updateStencil;
    StencilCompiler::Build(
        updateStencil,
        "UpdateStencil",
        domain,
        StencilConfiguration<Real, BlockSize<32, 4> >(),
        pack_parameters(
            // Input fields
            Param<qfine_, cIn>(qfine),
            Param<qcoarseold_, cIn>(qcoarseold),
            Param<qcoarsenew_, cIn>(qcoarsenew),
            // Output fields
            Param<q_, cInOut>(q)
        ),
        define_loops(
            define_sweep<cKIncrement>(
                define_stages(
                    StencilStage<UpdateStage, IJRange<cComplete,0,0,0,0>, KRange<FullDomain,0,0> >()
                )
            )
        )
    );

    // Convection object
    logfile << " - Initializing convection with nu = " << conf.nu() << std::endl;
    Convection convection(conf.gridSize(), conf.gridSize(), conf.gridSize(), conf.dx(), conf.nu(), conf.cx(), conf.cy(), conf.cz());

    // Fill initial solution
    SynchronizeHost(qinitial);
    fillQ(qinitial, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.); 
    SynchronizeDevice(qinitial);

    // Measure time required by convection
    double tauF = MPI_Wtime();
    convection.DoRK4(qinitial, qinitial, conf.dtFine(), conf.timeStepsFinePerTimeSlice());
    SynchronizeCUDA();
    tauF = MPI_Wtime() - tauF;
    double tauG = MPI_Wtime();
    convection.DoEuler(qinitial, qinitial, conf.dtCoarse(), conf.timeStepsCoarsePerTimeSlice());
    SynchronizeCUDA();
    tauG = MPI_Wtime() - tauG;

    const double tauRatio = tauG / tauF;
    const double Nit_Np = static_cast<double>(conf.kmax()) / commsize;
    const double maxSpeedup = 1. / (tauRatio * (1. + Nit_Np) + Nit_Np);

    // Fill again initial solution
    SynchronizeHost(qinitial);
    fillQ(qinitial, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.); 
    SynchronizeDevice(qinitial);

    // Prepare MPI objects
    MPI_Request reqSend, reqRecv;
    MPI_Status status;
    int mpiret;
    double * const pRecv = qinitial.storage().pStorageBase();
    double * const pSend = q.storage().pStorageBase();

    const IJKSize psize = q.storage().paddedSize();
    const int dataSize = psize.iSize() * psize.jSize() * psize.kSize();

    // Timers
    InternalTimer internalTimer;
    std::vector<double> times(6*conf.kmax() + 2);


    /**********************************
     * PARAREAL ALGORITHM STARTS HERE *
     **********************************/

    // Compute reference solution
    double eserial = MPI_Wtime();
    convection.DoRK4(qinitial, qreference, conf.dtFine(), (commrank+1)*conf.timeStepsFinePerTimeSlice());
    SynchronizeCUDA();
    eserial = MPI_Wtime() - eserial;
    matfile.addField(qreference);
    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize timing
    times[0] = MPI_Wtime();

    /*
     * Part 1: initialization
     */
    convection.DoEuler(qinitial, qinitial, conf.dtCoarse(), conf.timeStepsCoarsePerTimeSlice()*commrank);
    convection.DoEuler(qinitial, qcoarseold, conf.dtCoarse(), conf.timeStepsCoarsePerTimeSlice());

    convection.DoRK4(qinitial, qfine, conf.dtFine(), 1);

#ifdef ENABLE_INTERNAL_TIMING
    times[1] = MPI_Wtime();
#endif

    if (conf.mat())
    {
        matfile.addField(qinitial);
        matfile.addField(qcoarseold);
    }


    // Begin iteration
    for (int k = 0; k < conf.kmax(); ++k)
    {
        /*
         * Part 2: fine integration
         */

        // First step of fine propagation
        convection.DoRK4(qinitial, qfine, conf.dtFine(), 1);

        // Receive data in qinitial
        if (!isRoot && conf.async())
            MPI_Irecv(pRecv, dataSize, MPI_DOUBLE, commrank-1, k, MPI_COMM_WORLD, &reqRecv);

        // Remaining steps of fine propagation
        convection.DoRK4(qfine, qfine, conf.dtFine(), conf.timeStepsFinePerTimeSlice()-1);

        // Serialize
        if (conf.mat())
            matfile.addField(qfine, k);

#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        times[6*k + 2] = MPI_Wtime();
#endif

        /*
         * Part 3: coarse propagation and update
         */

        // Receive data
        if (!isRoot && conf.async())
            MPI_Wait(&reqRecv, &status);
        else if (!isRoot)
        {
            SynchronizeCUDA();
            MPI_Recv(pRecv, dataSize, MPI_DOUBLE, commrank-1, k, MPI_COMM_WORLD, &status);
            SynchronizeCUDA();
        }

        if (conf.mat())
            matfile.addField(qinitial, k);

#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        times[6*k + 3] = MPI_Wtime();
#endif

        // Coarse propagation
        convection.DoEuler(qinitial, qcoarsenew, conf.dtCoarse(), conf.timeStepsCoarsePerTimeSlice());

        if (conf.mat())
            matfile.addField(qcoarsenew, k);

#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        times[6*k + 4] = MPI_Wtime();
#endif

        // Wait for send
        if (!isLast && k > 0 && conf.async())
            MPI_Wait(&reqSend, &status);

#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        times[6*k + 5] = MPI_Wtime();
#endif

        // Update solution
        updateStencil.Apply();

#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        times[6*k + 6] = MPI_Wtime();
#endif

        if (conf.mat())
            matfile.addField(q, k);

        // Send solution to next process
        if (!isLast && conf.async())
            MPI_Isend(pSend, dataSize, MPI_DOUBLE, commrank+1, k, MPI_COMM_WORLD, &reqSend);
        else if (!isLast)
        {
            SynchronizeCUDA();
            MPI_Send(pSend, dataSize, MPI_DOUBLE, commrank+1, k, MPI_COMM_WORLD);
            SynchronizeCUDA();
        }

        // Swap qcoarse solutions
        qcoarsenew.SwapWith(qcoarseold);


        /*
         * Part 4: compute error
         */
#ifdef ENABLE_ERROR
        logfile << " - Error after iteration " << k << ": " << computeErrorReference(q, qreference)
                << std::endl;
#endif

        logfile << std::endl;

#ifdef ENABLE_INTERNAL_TIMING
        SynchronizeCUDA();
        times[6*k + 7] = MPI_Wtime();
#endif
    }

    // Get total timing
    double eparareal = MPI_Wtime() - times[0];

    // Wait for last send
    if (!isLast && conf.async())
        MPI_Wait(&reqSend, &status);


    /********************************
     * PARAREAL ALGORITHM ENDS HERE *
     ********************************/

    // Collect times
    for (int i = 6*conf.kmax()+1; i > 0; --i)
        times[i] = times[i] - times[i-1];

    logfile << "\n\nTimings:\n";
    for (int i = 0; i < 6*conf.kmax()+2; ++i)
        logfile << " - " << times[i] << "\n";

    std::vector<double> timesGlobal(isRoot ? (6*conf.kmax()+2) * (commsize) : 0);
    MPI_Gather(&times[0], 6*conf.kmax()+2, MPI_DOUBLE, &timesGlobal[0], 6*conf.kmax()+2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef ENABLE_INTERNAL_TIMING
    // Print times
    if (isRoot)
    {
        const int b = 6*conf.kmax() + 2;

        std::cout <<     "'Initialization             ':   ";
        for (int p = 0; p < commsize; ++p)
            std::cout << std::setw(12) << timesGlobal[p*b + 1] << "  ";
        std::cout << "\n";

        for (int k = 0; k < conf.kmax(); ++k)
        {
            std::cout << "'Fine propagation  , step " << k << " ':   ";
            for (int p = 0; p < commsize; ++p)
                std::cout << std::setw(12) << timesGlobal[p*b + 6*k + 2] << "  ";
            std::cout << "\n";

            std::cout << "'Receive data      , step " << k << " ':   ";
            for (int p = 0; p < commsize; ++p)
                std::cout << std::setw(12) << timesGlobal[p*b + 6*k + 3] << "  ";
            std::cout << "\n";

            std::cout << "'Coarse propagation, step " << k << " ':   ";
            for (int p = 0; p < commsize; ++p)
                std::cout << std::setw(12) << timesGlobal[p*b + 6*k + 4] << "  ";
            std::cout << "\n";

            std::cout << "'Send data         , step " << k << " ':   ";
            for (int p = 0; p < commsize; ++p)
                std::cout << std::setw(12) << timesGlobal[p*b + 6*k + 5] << "  ";
            std::cout << "\n";

            std::cout << "'Update solution   , step " << k << " ':   ";
            for (int p = 0; p < commsize; ++p)
                std::cout << std::setw(12) << timesGlobal[p*b + 6*k + 6] << "  ";
            std::cout << "\n";
        }
    }
#endif

    // Serializing last result
    matfile.addField("solution", q);

    // Writing error
    MPI_Barrier(MPI_COMM_WORLD);
    if (isLast)
    {
        double e = computeErrorReference(q, qreference);
        std::cout << "\n"
                  << "Error at end: " << e << "\n"
                  << "Speedup: " << eserial / eparareal << "\n"
                  << "Maximal speedup: " << maxSpeedup << std::endl;
    }

    // Closing log file
    logfile << " - Exiting. Goodbye\n\n";
    logfile << std::endl;

    // Finalize
    MPI_Finalize();

    return 0;
}

