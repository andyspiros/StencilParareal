#include <iostream>
#include <iomanip>
#include "RuntimeConfiguration.h"
#include "ConvectionSolution.h"
#include "mpi.h"
#include "Convection.h"
#include "SharedInfrastructure.h"
#include "StencilFramework.h"
#include "MatFile.h"

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

    double error = 0.;

    for (int i = jstart; i < iend; ++i)
        for (int j = jstart; j < jend; ++j)
            for (int k = kstart; k < kend; ++k)
            {
                double e = std::abs(field(i,j,k) - reference(i,j,k));
                error = std::max(error, e);
            }

#ifdef __CUDA_BACKEND__
    if (synchronizeField)
        field.SynchronizeDeviceStorage();
    if (synchronizeReference)
        reference.SynchronizeDeviceStorage();
#endif

    return error;
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int commsize, commrank;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
    const bool isRoot = commrank == 0;
    const bool isLast = commrank == commsize - 1;
    int mpiret;

    RuntimeConfiguration conf(argc, argv);

    // Compute timesteps
    const double dt = conf.endTime() / commsize;
    double dtFine = conf.dtFine();
    double dtCoarse = conf.dtCoarse();
    int timestepsFine = static_cast<int>(dt / dtFine + .5);
    int timestepsCoarse = static_cast<int>(dt / dtCoarse + .5);
    dtFine = dt / timestepsFine;
    dtCoarse = dt / timestepsCoarse;
    double cflFine = dtFine / (conf.dx()*conf.dx());
    double cflCoarse = dtCoarse / (conf.dx()*conf.dx());

    // Compute my start and end time
    const double timeStart = dt * commrank;
    const double timeEnd = dt * (commrank + 1);

    if (isRoot)
    std::cout << "Running with:\n"
        << " - diffusion coefficient: " << conf.nu() << "\n"
        << " - advection velocity in x: " << conf.cx() << "\n"
        << " - advection velocity in y: " << conf.cy() << "\n"
        << " - advection velocity in z: " << conf.cz() << "\n"
        << " - spatial discretization step: " << conf.dx() << "\n"
        << " - endtime: " << conf.endTime() << "\n"
        << " - number of time slices: " << commsize << "\n"
        << " - time slice size: " << dt << "\n"
        << " - CFL fine: " << cflFine << "\n"
        << " - CFL coarse: " << cflCoarse << "\n"
        << " - timestep size fine: " << dtFine << "\n"
        << " - timestep size coarse: " << dtCoarse << "\n"
        << " - timestep fine propagator: " << timestepsFine << "\n"
        << " - timestep coarse propagator: " << timestepsCoarse << "\n"
        << " - parareal iterations: " << conf.kmax() << "\n"
        << std::endl;

    // Reduce timestepsFine by one
    --timestepsFine;

    // Calculation domain and boundaries
    IJKSize domain; domain.Init(conf.gridSize(), conf.gridSize(), conf.gridSize()); 
    KBoundary kboundary; kboundary.Init(-convectionBoundaryLines, convectionBoundaryLines);

    // Instantiate fields
    ConvectionField q, qfromprevious, qfine, qcoarseold, qcoarsenew;
    ConvectionField qreference;
    q.Init("q", domain, kboundary);
    qfromprevious.Init("qfromprevious", domain, kboundary);
    qfine.Init("qfine", domain, kboundary);
    qcoarseold.Init("qcoarseold", domain, kboundary);
    qcoarsenew.Init("qcoarsenew", domain, kboundary);
    qreference.Init("qreference", domain, kboundary);

    // Create mat files
    std::ostringstream matfname;
    matfname << "Parareal_" << commrank << ".mat";
    MatFile matfile(matfname.str());

    // Create log file
    std::ostringstream logfname;
    logfname << "Parareal_" << commrank << ".log";
    std::ofstream logfile(logfname.str().c_str());

    logfile << "Rank " << commrank << " integrates from " << timeStart
        << " to " << timeEnd << "\n" << std::endl;

    // Field with initial conditions for root only,
    // otherwise, compute initial conditions in q
    ConvectionField* qinitial;
    if (isRoot)
    {
        qinitial = new ConvectionField;
        qinitial->Init("qinitial", domain, kboundary);
        SynchronizeHost(*qinitial);
        fillQ(*qinitial, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
        SynchronizeDevice(*qinitial);
        qfromprevious = *qinitial;
    }
    else
        fillQ(qfromprevious, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
    SynchronizeDevice(qfromprevious);

    SynchronizeCUDA();
    matfile.addField("qinitial", qfromprevious);
    logfile << " - Initial condition computed" << std::endl;

    // MPI requests for asynchronous communication
    MPI_Request reqRecv, reqSend;
    MPI_Status status;

    // Base storage for q and size of data
    const IJKSize& psize = q.storage().paddedSize();
    const int datasize = psize.iSize()*psize.jSize()*psize.kSize();
    double *qFromPreviousBuffer = qfromprevious.storage().pStorageBase();
    double *qToNextBuffer = q.storage().pStorageBase();
    logfile << " - qFromPreviousBuffer is " << qFromPreviousBuffer << std::endl;
    logfile << " - qToNextBuffer is " << qToNextBuffer << std::endl;

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

    // Convectionj object
    logfile << " - Initializing convection with nu = " << conf.nu() << std::endl;
    Convection convection(conf.gridSize(), conf.gridSize(), conf.gridSize(), conf.dx(), conf.nu(), conf.cx(), conf.cy(), conf.cz());

    // Timer
    InternalTimer internalTimer;

    // Fill reference solution
    SynchronizeDevice(qreference);
    logfile << " - Integrate reference to " << (commrank+1)*(timestepsFine+1)*dtFine << std::endl;
    convection.DoRK4(qfromprevious, qreference, dtFine, (commrank+1)*(timestepsFine+1));
    SynchronizeCUDA();
    SynchronizeHost(qreference);


    /**********************************
     * PARAREAL ALGORITHM STARTS HERE *
     **********************************/


    /*
     * 1. Initialization
     */
    internalTimer.start();
    convection.DoEuler(qfromprevious, qfromprevious, dtCoarse, timestepsCoarse*commrank);
    convection.DoEuler(qfromprevious, qcoarseold, dtCoarse, timestepsCoarse);
    SynchronizeCUDA();
    logfile << " - Initialization done in " << internalTimer.elapsed() << " msec" << std::endl;

    if (isRoot)
        std::cout << "Initialization done" << std::endl;

    logfile << "\n";

    /*
     * Begin of iteration
     */
    for (int k = 0; k < conf.kmax(); ++k)
    {
        logfile << "k = " << k << std::endl;
        
        if (isRoot)
            std::cout << "Begin of iteration " << k << std::endl;

        /*
         * 2. Fine propagation
         */
        internalTimer.start();
        convection.DoRK4Timestep(qfromprevious, qfine, dtFine);
        logfile << " - First timestep of fine propagator done in " << internalTimer.elapsed() << " msec" << std::endl;

        if (!isRoot)
        {
            internalTimer.start();
            mpiret = MPI_Irecv(qFromPreviousBuffer, datasize, MPI_DOUBLE, commrank-1, k, MPI_COMM_WORLD, &reqRecv);
            logfile << " - MPI_Irecv returned " << mpiret << " in " << internalTimer.elapsed() << " msec"<< std::endl;
        }

        internalTimer.start();
        convection.DoRK4(qfine, qfine, dtFine, timestepsFine);
        logfile << " - Fine propagation done in " << internalTimer.elapsed() << " msec" << std::endl;


        /**
         * 3. Coarse propagation
         */

        if (!isLast && k > 0)
        {
            internalTimer.start();
            mpiret = MPI_Wait(&reqSend, &status);
            logfile << " - MPI_Wait for send returned " << mpiret << " in " << internalTimer.elapsed() << " msec"<< std::endl;
        }

        if (!isRoot)
        {
            // Wait for receive to complete
            internalTimer.start();
            mpiret = MPI_Wait(&reqRecv, &status);
            logfile << " - MPI_Wait for recv returned " << mpiret << " in " << internalTimer.elapsed() << " msec"<< std::endl;

            // Integrate with coarse propagator
            internalTimer.start();
            convection.DoEuler(qfromprevious, qcoarsenew, dtCoarse, timestepsCoarse);
            logfile << " - Coarse propagator done in " << internalTimer.elapsed() << " msec" << std::endl;

            // Update solution
            internalTimer.start();
            updateStencil.Apply();
            logfile << " - Solution updated in " << internalTimer.elapsed() << " msec" << std::endl;

            // Store old coarse solution
            qcoarseold.SwapWith(qcoarsenew);
        }
        else
        {
            if(k == 0)
            {
                q.SwapWith(qfine);
                qToNextBuffer = q.storage().pStorageBase();
            }
        }

        // Send to next process
        if (!isLast)
        {
            internalTimer.start();
            mpiret = MPI_Isend(qToNextBuffer, datasize, MPI_DOUBLE, commrank+1, k, MPI_COMM_WORLD, &reqSend);
            logfile << " - MPI_Isend returned " << mpiret << " in " << internalTimer.elapsed() << " msec"<<  std::endl;
        }



        /**
         * 4. Error estimation
         */

#ifdef ENABLE_ERROR
        // Compute error
        SynchronizeCUDA();
        double error = computeErrorReference(q, qreference);
        matfile.addField("solution_at", q, k);
        matfile.addField("reference_at", qreference, k);
        SynchronizeDevice(q);

        logfile << " - Error at end of iteration " << k << ": " << error << "\n";
#endif
        logfile << std::endl;
    }

    logfile << "Parareal iteration finished. Cleaning up\n";

    // Wait for the last send
    if (!isLast)
    {
        internalTimer.start();
        mpiret = MPI_Wait(&reqSend, &status);
        logfile << " - MPI_Wait for last send returned " << mpiret << " in " << internalTimer.elapsed() << " msec"<< std::endl;
    }

    logfile << " - Saving solution and reference to file\n";
    matfile.addField("reference", qreference);
    matfile.addField("solution", q);

    // Close matfile
    matfile.close();

    // Delete root field
    if (isRoot)
        delete qinitial;

    logfile << " - Finalization\n";

    // Compute error
    if (isLast)
    {
        SynchronizeHost(q);
        double error = computeErrorReference(q, qreference);
        SynchronizeDevice(q);

        MPI_Send(&error, 1, MPI_DOUBLE, 0, conf.kmax(), MPI_COMM_WORLD);
    }
    if (isRoot)
    {
        double error;
        MPI_Status status;
        MPI_Recv(&error, 1, MPI_DOUBLE, commsize-1, conf.kmax(), MPI_COMM_WORLD, &status);

        std::cout << "Final error at time " << conf.endTime() << " is "
            << std::scientific << std::setprecision(10) << error << "\n";
    }

    // Finalize
    MPI_Finalize();

    logfile << " - Exiting. Goodbye\n\n";
    logfile << std::endl;

    return 0;
}

