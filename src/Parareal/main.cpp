#include <iostream>
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

class TTimer
{
public:

    TTimer()
    {
#ifdef ENABLE_TIMING
        std::cout << "Timing enabled\n";
#else
        std::cout << "Timing disabled\n";
#endif
    }

    void start()
    {
#ifdef ENABLE_TIMING
        cudaDeviceSynchronize();
        e_ = MPI_Wtime();
#endif
    }

    double elapsed()
    {
#ifdef ENABLE_TIMING
        cudaDeviceSynchronize();
        return (MPI_Wtime() - e_) * 1000.;
#else
        return 0.;
#endif
    }

private:
        double e_;

};

template<typename TDataField>
double computeErrorReference(const TDataField& field, const TDataField& reference)
{
    const IJKSize& size = field.calculationDomain();
    const int isize = size.iSize();
    const int jsize = size.jSize();
    const int ksize = size.kSize();

    double error = 0.;

    for (int i = 0; i < isize; ++i)
        for (int j = 0; j < jsize; ++j)
            for (int k = 0; k < ksize; ++k)
            {
                double e = std::abs(field(i,j,k) - reference(i,j,k));
                error = std::max(error, e);
            }
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
        qinitial->SynchronizeHostStorage();
        fillQ(*qinitial, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
        qinitial->SynchronizeDeviceStorage();
        qfromprevious = *qinitial;
    }
    else
        fillQ(qfromprevious, conf.nu(), conf.cx(), conf.cy(), conf.cz(), 0., 0., 1., 0., 1., 0., 1.);
    qfromprevious.SynchronizeDeviceStorage();

    cudaDeviceSynchronize();
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
    TTimer timer;

    // Fill reference solution
    qreference.SynchronizeDeviceStorage();
    logfile << " - Integrate reference to " << (commrank+1)*(timestepsFine+1)*dtFine << std::endl;
    convection.DoRK4(qfromprevious, qreference, dtFine, (commrank+1)*(timestepsFine+1));
    cudaDeviceSynchronize();
    qreference.SynchronizeHostStorage();


    /**********************************
     * PARAREAL ALGORITHM STARTS HERE *
     **********************************/


    /*
     * 1. Initialization
     */
    timer.start();
    convection.DoEuler(qfromprevious, qfromprevious, dtCoarse, timestepsCoarse*commrank);
    convection.DoEuler(qfromprevious, qcoarseold, dtCoarse, timestepsCoarse);
    cudaDeviceSynchronize();
    logfile << " - Initialization done in " << timer.elapsed() << " msec" << std::endl;

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
        timer.start();
        convection.DoRK4Timestep(qfromprevious, qfine, dtFine);
        logfile << " - First timestep of fine propagator done in " << timer.elapsed() << " msec" << std::endl;

        if (!isRoot)
        {
            timer.start();
            mpiret = MPI_Irecv(qFromPreviousBuffer, datasize, MPI_DOUBLE, commrank-1, k, MPI_COMM_WORLD, &reqRecv);
            logfile << " - MPI_Irecv returned " << mpiret << " in " << timer.elapsed() << " msec"<< std::endl;
        }

        timer.start();
        convection.DoRK4(qfine, qfine, dtFine, timestepsFine);
        logfile << " - Fine propagation done in " << timer.elapsed() << " msec" << std::endl;


        /**
         * 3. Coarse propagation
         */

        if (!isLast && k > 0)
        {
            timer.start();
            mpiret = MPI_Wait(&reqSend, &status);
            logfile << " - MPI_Wait for send returned " << mpiret << " in " << timer.elapsed() << " msec"<< std::endl;
        }

        if (!isRoot)
        {
            // Wait for receive to complete
            timer.start();
            mpiret = MPI_Wait(&reqRecv, &status);
            logfile << " - MPI_Wait for recv returned " << mpiret << " in " << timer.elapsed() << " msec"<< std::endl;

            // Integrate with coarse propagator
            timer.start();
            convection.DoEuler(qfromprevious, qcoarsenew, dtCoarse, timestepsCoarse);
            logfile << " - Coarse propagator done in " << timer.elapsed() << " msec" << std::endl;

            // Update solution
            timer.start();
            updateStencil.Apply();
            logfile << " - Solution updated in " << timer.elapsed() << " msec" << std::endl;

            // Store old coarse solution
            qcoarseold.SwapWith(qcoarsenew);
        }
        else
        {
            if(k == 0)
                q.SwapWith(qfine);
        }

        // Send to next process
        if (!isLast)
        {
            timer.start();
            mpiret = MPI_Isend(qToNextBuffer, datasize, MPI_DOUBLE, commrank+1, k, MPI_COMM_WORLD, &reqSend);
            logfile << " - MPI_Isend returned " << mpiret << " in " << timer.elapsed() << " msec"<<  std::endl;
        }



        /**
         * 4. Error estimation
         */

#ifdef ENABLE_ERROR
        // Compute error
        q.SynchronizeHostStorage();
        double error = computeErrorReference(q, qreference);
        q.SynchronizeDeviceStorage();

        logfile << " - Error at end of iteration " << k << ": " << error << "\n";
#endif
        logfile << std::endl;
    }

    logfile << "Parareal iteration finished. Cleaning up\n";

    // Wait for the last send
    if (!isLast)
    {
        timer.start();
        mpiret = MPI_Wait(&reqSend, &status);
        logfile << " - MPI_Wait for last send returned " << mpiret << " in " << timer.elapsed() << " msec"<< std::endl;
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

    // Finalize
    MPI_Finalize();

    logfile << " - Exiting. Goodbye\n\n";
    logfile << std::endl;

    return 0;
}

